# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Bernard UGUEN		: buguen@univ-rennes1.fr
#Mohamed LAARAIEDH	: mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import os
from numpy import *
from scipy import *
from scipy import optimize
from numpy.linalg import *
import cvxmod as cvxm
import cvxopt as cvxo
from CRBLocation import *

class ToALocation(object):
        """
        A ToALocation contains:
        1- a set of RadioNodes (RN) with associated position accuracies (RNQoS),
        2- a set of ToAs measurements (ToA) with associated STD (ToAStd) 
        
        This class manages the ToA-based localization techniques.
        MEMBERS:

                RN      : An Array that defines the Radio nodes implied in localization (coordiantes in meters)
                        : shape(RN)= (2 or 3,RNnum)
                RNQoS   : An Array that defines the precision of positions of RN (std in meters)
                        : shape(RNQoS)= (2 or 3, RNnum)
                ToA     : A measurement vector of ToA associated to RN (ToA values in seconds)
                        : shape(ToA)= (RNnum,1)
                ToAStd  : Associated STD of ToA (std in seconds)
                        : shape(ToAStd)= (RNnum,1)



        Provided Methods:
        
                info()                                          : Display information about the scenario

                LSToALocate(RN, ToA)                            : Applies Least Square approximation and estimate position
                WLSToALocate(RN, ToA, ToAStd)                   : Applies Weighted Least Square approximation and estimate position
                TSToALocation(P0, RN, ToA, ToAStd, Niter)       : Applies Taylor Series method and estimate position after Niter iterations
                AMLToALocation(P0, RN, ToA, ToAStd, Niter)      : This applies Aproximate Maximum Likelihood and estimate position after Niter
                TSMLToALocation(RN, ToA, ToAStd, Niter)         : This applies Two Stages Maximum Likelihood Method and estimates position
                
                ToAOptimizer(RN, ToA, ToAStd)                   : Defines the function to be optimized
                MLToALocate(P0, RN, ToA, ToAStd)                : Optimize ToAOptimizer() and estimate Position (P0:initial guess)

                CRBToALocate(self, P, RN, ToA, ToAStd)          : Compute the CRB in P for the given scenario
        """
	"""
        def __init__(self,RN, ToA, ToAStd):
                self.RN         = RN
                self.ToA        = ToA
                self.ToAStd     = ToAStd
	"""
        def __init__(self,RN):
                self.RN         = RN
	

        def info(self):
                """
                Dispaly scenario information
                """
                print "Reference Radio Nodes:\n", self.RN
                print "Measured ToA:\n", self.ToA
                print "STD of Measured ToA:\n", self.ToAStd

        def LSToALocate(self,RN, ToA):
                """
                This applies LS approximation on ToA to get position P.
                Return P
                """
                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                # Construct the vector K (see theory)
                RN2     = (sum(RN*RN,axis=0)).reshape(RNnum,1)
                k1      = RN2[1:RNnum,:]-RN2[0,0]                                       # first half of K
                
                RoA     = c*ToA                                         # Range of arrival (meters)
                RoA2    = (RoA*RoA).reshape(RNnum,1)
                k2      = RoA2[0,0]-RoA2[1:RNnum,:]                     # second half of K
                
                K       = k1+k2
                
                # Construct the matrix A (see theory)
                A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])

                # Apply LS operator
                P       = 0.5*dot(linalg.inv(dot(A.T,A)),dot(A.T,K))
                
                # Return the estimated position
                return P

        def  TLSToALocate(self, RN, ToA, ToAStd):
                """
                This applies TLS approximation on ToA to get position P.
                Return P
                """
        
                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                # Construct the vector K (see theory)
                RN2     = (sum(RN*RN,axis=0)).reshape(RNnum,1)
                k1      = RN2[1:RNnum,:]-RN2[0,0]                                       # first half of K
                
                RoA     = c*ToA                                         # Range of arrival (meters)
                RoA2    = (RoA*RoA).reshape(RNnum,1)
                k2      = RoA2[0,0]-RoA2[1:RNnum,:]                     # second half of K
                
                K       = k1+k2
                
                # Construct the matrix A (see theory)
                A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])

                A2   = dot(transpose(A),A)

                [U,S,V]=svd(A2)
                J = 1/S
                rA=rank(A)
                m,n=shape(A)
                f=0
                if log10(cond(A))>=c*max(ToAStd):
                        f=f+1
                        for i in range(n-rA):
                                u = where(J==max(J))
                                J[u] = 0

                A2i = dot(dot(V.T,diag(J)),U.T)
     
                P       = 0.5*dot(A2i,dot(A.T,K))

                return P

        def WLSToALocate(self, RN, ToA, ToAStd):
                """
                This applies WLS approximation on ToA assuming ToAStd to get position P.
                Return P
                """
                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                # Construct the vector K (see theory)
                RN2     = (sum(RN*RN,axis=0)).reshape(RNnum,1)
                k1      = RN2[1:RNnum,:]-RN2[0,0]                                       # first half of K
                
                RoA     = c*ToA                                         # Range of arrival (meters)
                RoAStd  = c*ToAStd
                RoA2    = (RoA*RoA).reshape(RNnum,1)
                k2      = RoA2[0,0]-RoA2[1:RNnum,:]                     # second half of K
                
                K       = k1+k2
                
                # Construct the matrix A (see theory)
                A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])
                
                # Construct the Covariance Matrix
                C       = diag((RoAStd[1:RNnum,0])**2)          
                
                # Apply LS operator
                P       =  0.5*dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),K))
        
                # Return the estimated position
                return P

        def  TWLSToALocate(self, RN, ToA, ToAStd):
                """
                This applies TWLS approximation on ToA to get position P.
                Return P
                """
        
                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                # Construct the vector K (see theory)
                RN2     = (sum(RN*RN,axis=0)).reshape(RNnum,1)
                k1      = RN2[1:RNnum,:]-RN2[0,0]                                       # first half of K
                
                RoA     = c*ToA                                         # Range of arrival (meters)
                RoAStd  = c*ToAStd
                RoA2    = (RoA*RoA).reshape(RNnum,1)
                k2      = RoA2[0,0]-RoA2[1:RNnum,:]                     # second half of K
                
                K       = k1+k2
                
                # Construct the matrix A (see theory)
                A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])
                
                # Construct the Covariance Matrix
                C       = diag((RoAStd[1:RNnum,0])**2)          

                A2   = dot(A.T,dot(linalg.inv(C),A))

                [U,S,V]=svd(A2)
                J = 1/S
                rA=rank(A)
                m,n=shape(A)
                f=0
                if log10(cond(A))>=c*max(ToAStd):
                        f=f+1
                        for i in range(n-rA):
                                u = where(J==max(J))
                                J[u] = 0

                A2i = dot(dot(V.T,diag(J)),U.T)
     
                P       = 0.5*dot(A2i,dot(dot(A.T,linalg.inv(C)),K))

                return P

        def TSToALocation(self, P0, RN, ToA, ToAStd, Niter):
                '''
                Applies Taylor Series method and estimate position after Niter iterations
                '''

                P       = P0                                            # Initialisation of P as equal to intial guess P0
                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                RoA     = c*ToA
                RoAStd  = c*ToAStd
                
                for i in arange(Niter):

                        # Construct the matrix A (see theory)
                        A       = ((outer(P,ones(RNnum))- RN)/sqrt(sum((outer(P,ones(RNnum))- RN)**2,axis=0))).T
                        
                        # Construct the Covariance Matrix
                        C       = diag((RoAStd[:,0])**2)

                        # Construct the vector D (see theory)
                        D       = RoA-sqrt((sum((outer(P,ones(RNnum))- RN)**2,axis=0)).reshape(shape(RoA)))

                        # construct the vector Delta (see theory)
                        Delta   = dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),D))

                        # update P
                        P       = P+Delta

                # Return the estimated position
                return P

        
        def ToAOptimizer1(self, P, RN, ToA, ToAStd):
                """
                This defines the ML function to be minimized
                """

                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                RoA     = c*ToA
                RoAStd  = c*ToAStd

                # construct the ML function to be minimized     
                RNmP    = RN - outer(P,ones(RNnum))
                mRNmP   = (sqrt(diag(dot(RNmP.T,RNmP)))).reshape(RNnum,1)
                tk      = (RoA-mRNmP)/mRNmP
                uk      = (tk * RNmP.T)/(RoAStd**2)
                msuk     = uk.sum(axis=0)
                #msuk    = sqrt(dot(suk,suk.T))

                
                return(msuk)

        def ToAOptimizer(self, P, RN, ToA, ToAStd):
                """
                This defines the ML function to be minimized
                """

                shRN    = shape(RN)                                     # shape of RN
                RNnum   = shRN[1]                                       # Number of reference nodes
                c       = 3e08                                          # Speed of light
                RoA     = c*ToA
                RoAStd  = c*ToAStd

                # construct the ML function to be minimized     
                RNmP    = RN - outer(P,ones(RNnum))
                mRNmP   = (sqrt(diag(dot(RNmP.T,RNmP)))).reshape(RNnum,1)
                tk      = (RoA-mRNmP)**2
                uk      = tk/(2*RoAStd**2)+log(sqrt(2*pi)*RoAStd)
                suk     = uk.sum(axis=0)
                #print suk
                msuk    = suk#sqrt(dot(suk,suk.T))

                
                return(msuk)
                

        def MLToALocate(self, P, P0, RN, ToA, ToAStd):
                """ 
                Optimization Routine
                """

                P       = optimize.fmin(self.ToAOptimizer,P0,args=(RN,ToA,ToAStd),xtol=1e-10,ftol=1e-10)

                return P.reshape(shape(P0))

        '''def SDPToALocate(self, RN, ToA, ToAStd):
                """
                Apply SDP approximation and localization
                """
                RN = cvxm.matrix(RN)
                ToA = cvxm.matrix(ToA)
                c       = 3e08                                          # Speed of light
                RoA     = c*ToA
                mtoa,ntoa=cvxm.size(RN)
                Im = cvxm.eye(mtoa)
                Y=cvxm.optvar('Y',mtoa+1,mtoa+1)
                t=cvxm.optvar('t',ntoa,1)
                prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
                prob.constr.append(Y>=0)
                prob.constr.append(Y[mtoa,mtoa]==1)
                for i in range(ntoa):
                    X0=cvxm.matrix([[Im, -cvxm.transpose(RN[:,i])],[-RN[:,i], cvxm.transpose(RN[:,i])*RN[:,i]]])
                    prob.constr.append(-t[i]<cvxm.trace(X0*Y)-RoA[i]**2)
                    prob.constr.append(t[i]>cvxm.trace(X0*Y)-RoA[i]**2)
                prob.solve()
                Pval=Y.value
                X_cvx=Pval[:2,-1]

                return X_cvx'''

        def SDPToALocate(self, RN, ToA, ToAStd):
                """
                Apply SDP approximation and localization
                """
                RN = cvxm.matrix(RN)
                ToA = cvxm.matrix(ToA)
                
                c       = 3e08                                          # Speed of light
                RoA     = c*ToA
                RoAStd  = c*ToAStd
                RoAStd = cvxm.matrix(RoAStd)
                mtoa,ntoa=cvxm.size(RN)
                Im = cvxm.eye(mtoa)
                Y=cvxm.optvar('Y',mtoa+1,mtoa+1)
                t=cvxm.optvar('t',ntoa,1)
                prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
                prob.constr.append(Y>=0)
                prob.constr.append(Y[mtoa,mtoa]==1)
                for i in range(ntoa):
                    X0=cvxm.matrix([[Im, -cvxm.transpose(RN[:,i])],[-RN[:,i], cvxm.transpose(RN[:,i])*RN[:,i]]])
                    prob.constr.append(-t[i]<(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                    prob.constr.append(t[i]>(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                prob.solve()
                Pval=Y.value
                X_cvx=Pval[:2,-1]

                return X_cvx
        
        def CRBToALocate(self, P, RN, ToAStd):
                """
                Compute the CRB in P for the given scenario
                """
                
                crlb=CRBLocation(RN)
                CRB=crlb.CRB_TOA_fim(P, RN, ToAStd)
                
                return sqrt(abs(CRB))



