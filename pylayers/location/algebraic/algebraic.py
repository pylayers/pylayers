# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#PYLAYERS is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#PYLAYERS is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with PYLAYERS.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Bernard UGUEN          : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################

import numpy as np
import scipy as sp
from scipy import optimize
import numpy.linalg as nplg
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from pylayers.util.geomutil import dist
import string

class algloc(object):
    """
    This class regroups all the algebraic localization sceanrios and
    techniques
    Attributes
    ----------
        nodes : dictionnary
        ldp : dictionnary
        c : speed of light
          float
        
    """

    def __init__(self, nodes, ldp):
        self.nodes = nodes
        self.ldp = ldp
        self.c = 0.3

    def info(self):
        """
        Dispaly scenario information
        """
        print "Speed of light : ", self.c
        print "Nodes : ", self.nodes
        print "Location dependant parameters : ", self.ldp

    def plot(self, rss, toa, tdoa):
        """ Plot scenario

        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean
        """

        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            fig = plt.figure()
            ax = Axes3D(fig)
            BN = self.nodes['BN']
            try:
                ax.plot(BN[0,:],BN[1,:],BN[2,:],'r*',zdir='z',\
                        label='blind node')
            except:
                plt.plot(BN[0,:],BN[1,:],'r*', label='blind node')
            if rss!=0:
                RN_RSS = self.nodes['RN_RSS']
                try:
                    ax.plot(RN_RSS[0,:],RN_RSS[1,:],RN_RSS[2,:],'ro',\
                            zdir='z', label='RSS node')
                except:
                    plt.plot(RN_RSS[0,:],RN_RSS[1,:],'ro',\
                            label='RSS node')
            if toa!=0:
                RN_TOA = self.nodes['RN_TOA']
                try:
                    ax.plot(RN_TOA[0,:],RN_TOA[1,:],RN_TOA[2,:],'gs',\
                            zdir='z', label='TOA node')
                except:
                    plt.plot(RN_TOA[0,:],RN_TOA[1,:], 'gs', \
                            label='TOA node')

            if toa!=0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                try:
                    ax.plot(RN_TDOA[0,:],RN_TDOA[1,:],RN_TDOA[2,:],'bD',\
                            zdir='z', label='TDOA node')
                    ax.plot(RNr_TDOA[0,:],RNr_TDOA[1,:],RNr_TDOA[2,:],\
                            'kD',zdir='z', label='Ref TDOA node')
                except:
                    plt.plot(RN_TDOA[0,:],RN_TDOA[1,:], 'bD', \
                            label='TDOA node')
                    plt.plot(RNr_TDOA[0,:],RNr_TDOA[1,:], 'kD', \
                            label='Ref TDOA node')

    def show(self, rss, toa, tdoa):
        """ Plot scenario

        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> S.show(1,1,1)
            >>> plt.show()
        """
        self.plot(rss, toa, tdoa)
        plt.legend(numpoints=1)
        plt.show()
        

    def get_range(self, Rest='mode'):
        """ Compute the range from RSS using the "Rest" estimator

        Parameters
        ----------
            Rest : string
        Returns
        -------
            Range : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> r_mode = S.get_range('mode')
            >>> r_median = S.get_range('median')
            >>> r_mean = S.get_range('mean')
        """
        RSS = self.ldp['RSS']
        RSS_std = self.ldp['RSS_std']
        RSS_np = self.ldp['RSS_np']
        d0 = self.ldp['d0']
        PL0 = self.ldp['PL0']
        S = (np.log(10)/10)*RSS_std/RSS_np
        M = (np.log(10)/10)*(PL0-RSS)/RSS_np+np.log(d0)
        if string.lower(Rest) == 'mode':
            Range = np.exp(M-S**2)
        elif string.lower(Rest) == 'median':
            Range = np.exp(M)
        elif string.lower(Rest) == 'mean':
            Range = np.exp(M+0.5*S**2)
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return Range

    def get_range_std(self, Rest='mode'):
        """
        Compute the RSS range standard deviation using the "Rest" \
        estimator
        Parameters
        ----------
            Rest : string
        Returns
        -------
            Range_std : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> rs_mode = S.get_range_std('mode')
            >>> rs_median = S.get_range_std('median')
            >>> rs_mean = S.get_range_std('mean')
        """
        RSS = self.ldp['RSS']
        RSS_std = self.ldp['RSS_std']
        RSS_np = self.ldp['RSS_np']
        d0 = self.ldp['d0']
        PL0 = self.ldp['PL0']
        S = (np.log(10)/10)*RSS_std/RSS_np
        M = (np.log(10)/10)*(PL0-RSS)/RSS_np+np.log(d0)
        if string.lower(Rest) == 'mode':
            Range_std = np.sqrt((np.exp(2*M-2*S**2))*(1-np.exp(-S**2)))
        elif string.lower(Rest) == 'median':
            Range_std = np.sqrt((np.exp(2*M+S**2))*(np.exp(S**2)-1))
        elif string.lower(Rest) == 'mean':
            Range_std = np.sqrt((np.exp(2*M+3*S**2))*(np.exp(S**2)-1))
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return Range_std
        
    def ls_locate(self, rss, toa, tdoa, Rest):
        """
        This method applies least squares (LS) approximation to get
        position P.
        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.ls_locate(1, 0, 0, 'mode')
            >>> P_toa = S.ls_locate(0, 1, 0, 'mode')
            >>> P_tdoa = S.ls_locate(0, 0, 1, 'mode')
            >>> P_rsstoa = S.ls_locate(1, 1, 0, 'mode')
            >>> P_rsstdoa = S.ls_locate(1, 0, 1, 'mode')
            >>> P_toatdoa = S.ls_locate(0, 1, 1, 'mode')
            >>> P_rsstoatdoa = S.ls_locate(1, 1, 1, 'mode')
            
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K = k1-k2
                # Construct the matrix A (see theory)
                A = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                # Construct the vector K (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.c * TOA
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = P.reshape(np.shape(RN_TOA[:,0:1]))
                
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                # Construct the vector K (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.get_range(Rest)
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TOA.reshape(np.shape(K_RSS)[0],1)))
                A = np.vstack((A_RSS, A_TOA))
                P = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif rss==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2 = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1 + k2
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1_TOA = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2_TOA = RDOA2
                K_TDOA = k1_TOA-k2_TOA
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_TOA.reshape(np.shape(K_TOA)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_TOA, A_TDOA))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)

            elif toa == 0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_RSS, A_TDOA))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)

            else:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2 = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1 + k2
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1_TOA = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2_TOA = RDOA2
                K_TDOA = k1_TOA-k2_TOA
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((np.vstack((K_RSS.reshape(np.shape(K_RSS)\
                    [0],1), K_TOA.reshape(np.shape(K_TOA)[0],1))),\
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((np.vstack((A_RSS, A_TOA)),A_TDOA))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            return P


    def wls_locate(self, rss, toa, tdoa, Rest):
        """
        This method applies weighted least squares (WLS) approximation
        to get position P.
        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.wls_locate(1, 0, 0, 'mode')
            >>> P_toa = S.wls_locate(0, 1, 0, 'mode')
            >>> P_tdoa = S.wls_locate(0, 0, 1, 'mode')
            >>> P_rsstoa = S.wls_locate(1, 1, 0, 'mode')
            >>> P_rsstdoa = S.wls_locate(1, 0, 1, 'mode')
            >>> P_toatdoa = S.wls_locate(0, 1, 1, 'mode')
            >>> P_rsstoatdoa = S.wls_locate(1, 1, 1, 'mode')
            
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K = k1-k2
                # Construct the matrix A (see theory)
                A = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                RDOA.reshape(np.shape(TDOA)[0],1)))
                # Construct the Covariance Matrix
                C = np.diag(RDOA_std[:] ** 2)
                # Apply LS operator
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C)\
                ,A))),np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                # Construct the vector K (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.c * TOA
                ROA_std = self.c * TOA_std
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Construct the Covariance Matrix
                C = np.diag(ROA_std[1:] ** 2)
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),\
                A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = P.reshape(np.shape(RN_TOA[:,0:1]))
                
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                # Construct the vector K (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.get_range(Rest)
                ROA_std = self.get_range_std(Rest)
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Construct the Covariance Matrix
                C = np.diag((ROA_std[1:]) ** 2)
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),\
                A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA_RSS_std = self.get_range_std(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA_TOA_std = self.c * TOA_std
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TOA.reshape(np.shape(K_RSS)[0],1)))
                A = np.vstack((A_RSS, A_TOA))
                C = np.diag(np.hstack((ROA_RSS_std[1:]**2,\
                    ROA_TOA_std[1:]**2)))
                P = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),\
                A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif rss==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA_TOA_std = self.c * TOA_std
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_TOA.reshape(np.shape(K_TOA)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_TOA, A_TDOA))
                C = np.diag(np.hstack((ROA_TOA_std[1:]**2,\
                    RDOA_std[:]**2)))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,np.dot(nplg.inv(C),\
                    A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)

            elif toa == 0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA_RSS_std = self.get_range_std(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_RSS, A_TDOA))
                C = np.diag(np.hstack((ROA_RSS_std[1:]**2,\
                    RDOA_std[:]**2)))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,np.dot(nplg.inv(C),\
                    A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)

            else:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA_RSS_std = self.get_range_std(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA_TOA_std = self.c * TOA_std
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((np.vstack((K_RSS.reshape(np.shape(K_RSS)\
                    [0],1), K_TOA.reshape(np.shape(K_TOA)[0],1))),\
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((np.vstack((A_RSS, A_TOA)),A_TDOA))
                C = np.diag(np.hstack((np.hstack((ROA_RSS_std[1:]**2,\
                    ROA_TOA_std[1:]**2)),RDOA_std[:]**2)))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,np.dot(nplg.inv(C),\
                    A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            return P

    def ml_function(self, P, rss, toa, tdoa):
        """
        This defines the ML function to be minimized if ML estimator
        is used
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            ML : numpy.ndarray
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                RNnum = shRN[1]
                RDOA = self.c*TDOA
                RDOA_std = self.c*TDOA_std
                # construct the ML function to be minimized
                RN1mP = RN_TDOA - np.outer(P, np.ones(RNnum))
                mRN1mP = np.sqrt(np.diag(np.dot(RN1mP.T, RN1mP)))
                RN2mP = RNr_TDOA - np.outer(P, np.ones(RNnum))
                mRN2mP = np.sqrt(np.diag(np.dot(RN2mP.T, RN2mP)))
                rRDOA = mRN1mP - mRN2mP
                tk = (RDOA - rRDOA) ** 2 / (2 * RDOA_std ** 2)
                ML = tk.sum(axis=0)
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RNnum = np.shape(RN_TOA)[1]
                ROA = self.c * TOA
                ROA_std = self.c * TOA_std
                # construct the ML function to be minimized
                RNmP = RN_TOA - np.outer(P, np.ones(RNnum))
                mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
                tk = (ROA-mRNmP) ** 2
                uk = tk/(2*ROA_std**2)+np.log(np.sqrt(2*np.pi)*ROA_std)
                ML = uk.sum(axis=0)
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RNnum = np.shape(RN_RSS)[1]
                S = (np.log(10)/10)*RSS_std/RSS_np
                M = (np.log(10)/10)*(PL0-RSS)/RSS_np+np.log(d0)
                # construct the ML function to be minimized
                RNmP = RN_RSS-np.outer(P, np.ones(RNnum))
                mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
                tk = (M-S**2-np.log(mRNmP))**2
                uk = tk/(2*S**2)
                ML = uk.sum(axis=0)
            elif rss==0:
                ML = self.ml_function(P,0,1,0)+self.ml_function(P,0,0,1)
            elif toa==0:
                ML = self.ml_function(P,1,0,0)+self.ml_function(P,0,0,1)
            elif tdoa==0:
                ML = self.ml_function(P,1,0,0)+self.ml_function(P,0,1,0)
            else:
                ML =self.ml_function(P,1,0,0)+self.ml_function(P,0,1,0)\
                    + self.ml_function(P,0,0,1)

            return ML

    def ml_locate(self, P0, rss, toa, tdoa):
        """
        This applies ML estimator to compute position P
        Parameters
        ----------
            P0 : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.ml_locate(BN0, 1, 0, 0)
            >>> P_toa = S.ml_locate(BN0, 0, 1, 0)
            >>> P_tdoa = S.ml_locate(BN0, 0, 0, 1)
            >>> P_rsstoa = S.ml_locate(BN0, 1, 1, 0)
            >>> P_rsstdoa = S.ml_locate(BN0, 1, 0, 1)
            >>> P_toatdoa = S.ml_locate(BN0, 0, 1, 1)
            >>> P_rsstoatdoa = S.ml_locate(BN0, 1, 1, 1)
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            P = optimize.fmin(self.ml_function,P0,\
                args=(rss,toa,tdoa),disp=0)
            return P.reshape(np.shape(P0))

    def fim(self, P, rss, toa, tdoa):
        """
        Compute the fisher information matrix in P for the given
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            FIM : numpy.ndarray
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                RNnum = shRN[1]
                RDOA = self.c*TDOA
                RDOA_std = self.c*TDOA_std
                FIM = np.zeros((np.shape(RN_TDOA)[0],\
                        np.shape(RN_TDOA)[0]))
                for i in range(np.shape(RN_TDOA)[1]):
                    PmRN = np.sqrt(np.dot((P-RN_TDOA[:,i:i+1]).T,\
                        P-RN_TDOA[:,i:i+1]))
                    PmRNr = np.sqrt(np.dot((P-RNr_TDOA[:,i:i+1]).T,\
                        P-RNr_TDOA[:,i:i+1]))
                    FIM += (1/RDOA_std[i]**2)*\
                            np.dot((P-RN_TDOA[:,i:i+1])/PmRN-\
                            (P-RNr_TDOA[:,i:i+1])/PmRNr,\
                            ((P-RN_TDOA[:,i:i+1])/PmRN-\
                            (P-RNr_TDOA[:,i:i+1])/PmRNr).T)
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RNnum = np.shape(RN_TOA)[1]
                ROA = self.c * TOA
                ROA_std = self.c * TOA_std
                FIM = np.zeros((np.shape(RN_TOA)[0],np.shape(RN_TOA)[0]))
                for i in range(np.shape(RN_TOA)[1]):
                    FIM += np.dot((P-RN_TOA[:,i:i+1]),\
                        (P-RN_TOA[:,i:i+1]).T)/((ROA_std[i]**2)*\
                        np.dot((P-RN_TOA[:,i:i+1]).T,\
                        (P-RN_TOA[:,i:i+1])))
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RNnum = np.shape(RN_RSS)[1]
                S = (np.log(10)/10)*RSS_std/RSS_np
                FIM = np.zeros((np.shape(RN_RSS)[0],\
                    np.shape(RN_RSS)[0]))
                for i in range(np.shape(RN_RSS)[1]):
                    FIM += np.dot((P-RN_RSS[:,i:i+1]),\
                        (P-RN_RSS[:,i:i+1]).T)/((S[0]**2)*\
                        (np.dot((P-RN_RSS[:,i:i+1]).T,\
                        (P-RN_RSS[:,i:i+1])))**2)
            elif rss==0:
                FIM = self.fim(P,0,1,0)+self.fim(P,0,0,1)
            elif toa==0:
                FIM = self.fim(P,1,0,0)+self.fim(P,0,0,1)
            elif tdoa==0:
                FIM = self.fim(P,1,0,0)+self.fim(P,0,1,0)
            else:
                FIM = self.fim(P,1,0,0)+self.fim(P,0,1,0)+\
                    self.fim(P,0,0,1)

            return FIM

    def crb(self, P, rss, toa, tdoa):
        """
        This method compute the cramer rao bound (CRB) at position P.
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            CRB : float

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> crb_rss = S.crb(BN, 1, 0, 0)
            >>> crb_toa = S.crb(BN, 0, 1, 0)
            >>> crb_tdoa = S.crb(BN, 0, 0, 1)
            >>> crb_rsstoa = S.crb(BN, 1, 1, 0)
            >>> crb_rsstdoa = S.crb(BN, 1, 0, 1)
            >>> crb_toatdoa = S.crb(BN, 0, 1, 1)
            >>> crb_rsstoatdoa = S.crb(BN, 1, 1, 1)
        """

        FIM = self.fim(P, rss, toa, tdoa)
        CRB = np.sqrt(np.trace(nplg.inv(FIM)))
        return CRB

def scenario():
    """
    This method is not a class member, it defines a sample scenario
    Returns
    -------
        CRB : float

    Examples
    --------
    .. plot::
        :include-source:
        
        >>> nRN = 4
        >>> dim = 3 # 2 for 2D, 3 for 3D
        >>> L = 20.
        >>> c = 0.3
        >>> BN = L*sp.rand(dim,1)
        >>> BN0 = L*sp.rand(dim,1)
        >>> RN_TOA = L*sp.rand(dim,nRN)
        >>> RN_RSS = L*sp.rand(dim,nRN)
        >>> RN_TDOA = L*sp.rand(dim,nRN)
            
        >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
        >>> TOF = d_TOA/c # actual TOA
        >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
        >>> TOA = TOF + TOA_std

        >>> RSS_std = 0.001 * np.ones(nRN)
        >>> RSS_np = 2.645 * np.ones(nRN)
        >>> PL0 = 34.7*np.ones(nRN)
        >>> d0 = 1.
        >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
        >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
        >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
        >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
        >>> d = dist(RN_TDOA,BN,0)
        >>> dr = dist(RNr_TDOA,BN,0)
        >>> TDOF = (d-dr)/c # actual TDOA
        >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
        >>> TDOA = TDOF + TDOA_std

        >>> nodes={}
        >>> nodes['BN']= BN
        >>> nodes['RN_RSS']= RN_RSS
        >>> nodes['RN_TOA']= RN_TOA
        >>> nodes['RN_TDOA']= RN_TDOA
        >>> nodes['RNr_TDOA']= RNr_TDOA

        >>> ldp={}
        >>> ldp['RSS'] = RSS
        >>> ldp['RSS_std'] = RSS_std
        >>> ldp['RSS_np'] = RSS_np
        >>> ldp['d0'] = d0
        >>> ldp['PL0'] = PL0
        >>> ldp['TOA'] = TOA
        >>> ldp['TOA_std'] = TOA_std
        >>> ldp['TDOA'] = TDOA
        >>> ldp['TDOA_std'] = TDOA_std
            
        >>> print 'Nodes'
        >>> print nodes
        >>> print 'LDPs'
        >>> print ldp
        >>> print 'BN0:initial guess for ML estimator'
        >>> print BN0
        """
    nRN = 4
    dim = 3 # 2 for 2D, 3 for 3D
    L = 20.
    c = 0.3
    BN = L*sp.rand(dim,1)
    BN0 = L*sp.rand(dim,1)
    RN_TOA = L*sp.rand(dim,nRN)
    RN_RSS = L*sp.rand(dim,nRN)
    RN_TDOA = L*sp.rand(dim,nRN)
            
    d_TOA = dist(RN_TOA,BN,0) # actual distances
    TOF = d_TOA/c # actual TOA
    TOA_std = 0.001/c*np.ones(np.shape(TOF))
    TOA = TOF + TOA_std

    RSS_std = 0.001 * np.ones(nRN)
    RSS_np = 2.645 * np.ones(nRN)
    PL0 = 34.7*np.ones(nRN)
    d0 = 1.
    d_RSS = dist(RN_RSS,BN,0) # actual distances
    X = RSS_std * np.random.randn(np.shape(PL0)[0])
    RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
    RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
    d = dist(RN_TDOA,BN,0)
    dr = dist(RNr_TDOA,BN,0)
    TDOF = (d-dr)/c # actual TDOA
    TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
    TDOA = TDOF + TDOA_std

    nodes={}
    nodes['BN']= BN
    nodes['RN_RSS']= RN_RSS
    nodes['RN_TOA']= RN_TOA
    nodes['RN_TDOA']= RN_TDOA
    nodes['RNr_TDOA']= RNr_TDOA

    ldp={}
    ldp['RSS'] = RSS
    ldp['RSS_std'] = RSS_std
    ldp['RSS_np'] = RSS_np
    ldp['d0'] = d0
    ldp['PL0'] = PL0
    ldp['TOA'] = TOA
    ldp['TOA_std'] = TOA_std
    ldp['TDOA'] = TDOA
    ldp['TDOA_std'] = TDOA_std

    return nodes, ldp, BN0

        
    

if __name__ == "__main__":
    doctest.testmod()

        
                
                
  
