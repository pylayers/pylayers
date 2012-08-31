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
from numpy import *
from scipy import *
from scipy import optimize
from numpy.linalg import *
import cvxmod as cvxm
import cvxopt as cvxo
from string import *
from CRBLocation import *
from RSSLocation import *
from ToALocation import *
from TDoALocation import *

"""
                RN      : An Array that defines the Radio nodes implied in localization (coordiantes in meters)
                        : shape(RN)= (2 or 3,RNnum)
                RNQoS   : An Array that defines the precision of positions of RN (std in meters)
                        : shape(RNQoS)= (2 or 3, RNnum)
                ToA     : A measurement vector of ToA associated to RN (ToA values in seconds)
                        : shape(ToA)= (RNnum,1)
                ToAStd  : Associated STD of ToA (std in seconds)
                        : shape(ToAStd)= (RNnum,1)
"""
class HDFLocation(object):
        """
        This class regroups methods of localization using HDF of RSSI, TOA, and TDOA 
        """

        
        def __init__(self, RN_RSS, RN_ToA, RN_TDoA):
                self.RN_RSS         = RN_RSS
                self.RN_ToA         = RN_ToA
                self.RN_TDoA         = RN_TDoA
                

        def LSHDFLocate(self, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
                c       = 3e08
                RSSL=RSSLocation(RN_RSS)
                if RN_RSS==None:
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)
                    
                
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                
                        
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        # solution
                        K=vstack((K_ToA, K_TDoA))
                        A=vstack((A_ToA, A_TDoA))
                        P       = 0.5*dot(linalg.inv(dot(A.T,A)),dot(A.T,K))
                        return P[:2,:]

                elif RN_ToA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)
                    
                
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                
                        
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        # solution
                        K=vstack((K_RSS, K_TDoA))
                        A=vstack((A_RSS, A_TDoA))
                        P       = 0.5*dot(linalg.inv(dot(A.T,A)),dot(A.T,K))
                        return P[:2,:]
                
                elif RN_TDoA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0])
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0])
                        # solution
                        K=vstack((K_RSS, K_ToA))
                        A=vstack((A_RSS, A_ToA))
                        P       = 0.5*dot(linalg.inv(dot(A.T,A)),dot(A.T,K))
                        return P

                else:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA         = c*ToA                                  
                        RoA2        = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))

                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_TDoA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)
                    
                
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                
                        
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        # solution
                        K=vstack((vstack((K_RSS,K_ToA)), K_TDoA))
                        A=vstack((vstack((A_RSS,A_ToA)), A_TDoA))
                        P       = 0.5*dot(linalg.inv(dot(A.T,A)),dot(A.T,K))
                        return P[:2,:]

        def TLSHDFLocate(self, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
                c       = 3e08
                RSSL=RSSLocation(RN_RSS)
                if RN_RSS==None:
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)
                    
                
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                
                        
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        # solution
                        K=vstack((K_ToA, K_TDoA))
                        A=vstack((A_ToA, A_TDoA))
                        A2   = dot(transpose(A),A)

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(A.T,K))
                        
                        return P[:2,:]

                elif RN_ToA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)
                    
                
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                
                        
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        # solution
                        K=vstack((K_RSS, K_TDoA))
                        A=vstack((A_RSS, A_TDoA))
                        A2   = dot(transpose(A),A)

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(A.T,K))
                        return P[:2,:]
                
                elif RN_TDoA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0])
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0])
                        # solution
                        K=vstack((K_RSS, K_ToA))
                        A=vstack((A_RSS, A_ToA))
                        A2   = dot(transpose(A),A)

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(A.T,K))
                        return P

                else:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))

                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_TDoA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)
                    
                
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                
                        
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        # solution
                        K=vstack((vstack((K_RSS,K_ToA)), K_TDoA))
                        A=vstack((vstack((A_RSS,A_ToA)), A_TDoA))
                        A2   = dot(transpose(A),A)

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(A.T,K))
                        return P[:2,:]


        def WLSHDFLocate(self, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
                c       = 3e08
                RSSL=RSSLocation(RN_RSS)
                if RN_RSS==None:
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA
                        RoAStd     = c*ToAStd
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        C_ToA       = RoAStd[1:RNnum_ToA,0]**2
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoAStd    = c*TDoAStd
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)               
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        C_TDoA       = RDoAStd[:,0]**2
                        # solution
                        K=vstack((K_ToA, K_TDoA))
                        A=vstack((A_ToA, A_TDoA))
                        C=diag(hstack((C_ToA, C_TDoA)))
                        P       = 0.5*dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),K))
                        return P[:2,:]

                elif RN_ToA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoAStd      = RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        C_RSS       = RoAStd[1:RNnum_RSS,0]**2
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoAStd    = c*TDoAStd
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)               
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        C_TDoA       = RDoAStd[:,0]**2
                        # solution
                        K=vstack((K_RSS, K_TDoA))
                        A=vstack((A_RSS, A_TDoA))
                        C=diag(hstack((C_RSS, C_TDoA)))
                        P       = 0.5*dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),K))
                        return P[:2,:]
                
                elif RN_TDoA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoAStd      = RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        C_RSS       = RoAStd[1:RNnum_RSS,0]**2
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoAStd     = c*ToAStd
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        C_ToA       = RoAStd[1:RNnum_ToA,0]**2
                        # solution
                        K=vstack((K_RSS, K_ToA))
                        A=vstack((A_RSS, A_ToA))
                        C=diag(hstack((C_RSS, C_ToA)))
                        P       = 0.5*dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),K))
                        return P

                else:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoAStd      = RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        C_RSS       = RoAStd[1:RNnum_RSS,0]**2
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA         = c*ToA                                  
                        RoAStd     = c*ToAStd
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        C_ToA       = RoAStd[1:RNnum_ToA,0]**2

                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_TDoA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoAStd    = c*TDoAStd
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)               
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        C_TDoA       = RDoAStd[:,0]**2
                        # solution
                        K=vstack((vstack((K_RSS,K_ToA)), K_TDoA))
                        A=vstack((vstack((A_RSS,A_ToA)), A_TDoA))
                        C=diag(hstack((hstack((C_RSS,C_ToA)), C_TDoA)))
                        P       = 0.5*dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),K))
                        return P[:2,:]

        def TWLSHDFLocate(self, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
                c       = 3e08
		try:
	                RSSL=RSSLocation(RN_RSS)
		except :
			pass
                if RN_RSS==None:
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoAStd     = c*ToAStd
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        C_ToA       = RoAStd[1:RNnum_ToA,0]**2
                        
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoAStd    = c*TDoAStd
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)               
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        C_TDoA       = RDoAStd[:,0]**2
                        # solution
                        K=vstack((K_ToA, K_TDoA))
                        A=vstack((A_ToA, A_TDoA))
                        C=diag(hstack((C_ToA, C_TDoA)))
                        A2   = dot(A.T,dot(linalg.inv(C),A))

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                	
#                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
#                                f=f+1
#                                for i in range(n-rA):
#                                        u = where(J==max(J))
#                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(dot(A.T,linalg.inv(C)),K))
                        
                        return P[:2,:]

                elif RN_ToA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoAStd      = RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        C_RSS       = RoAStd[1:RNnum_RSS,0]**2
                        
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoAStd    = c*TDoAStd
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)               
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        C_TDoA       = RDoAStd[:,0]**2
                        # solution
                        K=vstack((K_RSS, K_TDoA))
                        A=vstack((A_RSS, A_TDoA))
                        C=diag(hstack((C_RSS, C_TDoA)))
                        A2   = dot(A.T,dot(linalg.inv(C),A))

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(dot(A.T,linalg.inv(C)),K))
                        return P[:2,:]
                
                elif RN_TDoA==None:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoAStd      = RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        C_RSS       = RoAStd[1:RNnum_RSS,0]**2
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoAStd     = c*ToAStd
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        C_ToA       = RoAStd[1:RNnum_ToA,0]**2
                        # solution
                        K=vstack((K_RSS, K_ToA))
                        A=vstack((A_RSS, A_ToA))
                        C=diag(hstack((C_RSS, C_ToA)))
                        A2   = dot(A.T,dot(linalg.inv(C),A))

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(dot(A.T,linalg.inv(C)),K))
                        return P

                else:
                        # for RSS
                        shRN_RSS    = shape(RN_RSS)                                     
                        RNnum_RSS   = shRN_RSS[1]                                    
                                                               
                        
                        RN_RSS2     = (sum(RN_RSS*RN_RSS,axis=0)).reshape(RNnum_RSS,1)      
                        RoA         = RSSL.getRange(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)        # RSS based Ranges (meters)
                        RoAStd      = RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
                        RoA2        = (RoA*RoA).reshape(RNnum_RSS,1)
                        K_RSS       = RN_RSS2[1:RNnum_RSS,:]-RN_RSS2[0,0] + RoA2[0,0]-RoA2[1:RNnum_RSS,:]
                
                        A_RSS       = hstack((RN_RSS[:,1:RNnum_RSS].T - RN_RSS[:,0].reshape(1,shRN_RSS[0]), zeros((RNnum_RSS-1,1))))
                        C_RSS       = RoAStd[1:RNnum_RSS,0]**2
                        
                        # for ToA
                        shRN_ToA    = shape(RN_ToA)                                     
                        RNnum_ToA   = shRN_ToA[1]                                    
                                                               
                        
                        RN_ToA2     = (sum(RN_ToA*RN_ToA,axis=0)).reshape(RNnum_ToA,1)      
                        RoA     = c*ToA                                  
                        RoAStd     = c*ToAStd
                        RoA2    = (RoA*RoA).reshape(RNnum_ToA,1)
                        K_ToA       = RN_ToA2[1:RNnum_ToA,:]-RN_ToA2[0,0] + RoA2[0,0]-RoA2[1:RNnum_ToA,:]
                
                        A_ToA       = hstack((RN_ToA[:,1:RNnum_ToA].T - RN_ToA[:,0].reshape(1,shRN_ToA[0]), zeros((RNnum_ToA-1,1))))
                        C_ToA       = RoAStd[1:RNnum_ToA,0]**2
                        # for TDoA
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_TDoA[:,0:1]*ones((1,RNnum_TDoA))
                        RDoA    = c*TDoA
                        RDoAStd    = c*TDoAStd
                        RDoA2   = (RDoA*RDoA).reshape(RNnum_TDoA,1)               
                        K_TDoA       = (sum((RN_TDoA-RN_TDoA2)*(RN_TDoA-RN_TDoA2),axis=0)).reshape(RNnum_TDoA,1)-RDoA2
                        A_TDoA       = hstack((RN_TDoA.T - RN_TDoA2.T,0.5*RoA[0,0]*RDoA))
                        C_TDoA       = RDoAStd[:,0]**2
                        # solution
                        K=vstack((vstack((K_RSS,K_ToA)), K_TDoA))
                        A=vstack((vstack((A_RSS,A_ToA)), A_TDoA))
                        C=diag(hstack((hstack((C_RSS,C_ToA)), C_TDoA)))
                        A2   = dot(A.T,dot(linalg.inv(C),A))

                        [U,S,V]=svd(A2)
                        J = 1/S
                        rA=rank(A)
                        m,n=shape(A)
                        f=0
                
                        if log10(cond(A))>=max(RSSL.getRangeStd(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
                                f=f+1
                                for i in range(n-rA):
                                        u = where(J==max(J))
                                        J[u] = 0

                        A2i = dot(dot(V.T,diag(J)),U.T)
     
                        P       = 0.5*dot(A2i,dot(dot(A.T,linalg.inv(C)),K))
                        return P[:2,:]

        def HDFOptimizer(self, P, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
               
                RSSL=RSSLocation(RN_RSS)
                TOAL=ToALocation(RN_ToA)
                TDOAL=TDoALocation(RN_TDoA)


                if RN_RSS==None:
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        fopt=TOAL.ToAOptimizer(P, RN_ToA, ToA, ToAStd) + TDOAL.TDoAOptimizer(P, RN_TDoA, RN_TDoA2, TDoA, TDoAStd)
                        return fopt       

                elif RN_ToA==None:
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        fopt=RSSL.DRSSOptimizer(P, RN_RSS, PL0, d0, RSS, RSSnp, RSSStd) + TDOAL.TDoAOptimizer(P, RN_TDoA, RN_TDoA2, TDoA, TDoAStd)
                        return fopt
                
                elif RN_TDoA==None:
                        fopt=RSSL.DRSSOptimizer(P, RN_RSS, PL0, d0, RSS, RSSnp, RSSStd) + TOAL.ToAOptimizer(P, RN_ToA, ToA, ToAStd)
                        return fopt

                else:
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_TDoA[:,0:1]*ones((1,RNnum_TDoA))
                        fopt=RSSL.DRSSOptimizer(P, RN_RSS, PL0, d0, RSS, RSSnp, RSSStd) + TOAL.ToAOptimizer(P, RN_ToA, ToA, ToAStd) + TDOAL.TDoAOptimizer(P, RN_TDoA, RN_TDoA2, TDoA, TDoAStd)
                        return fopt
                
        def MLHDFLocate(self, P, P0, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
                
                P       = optimize.fmin(self.HDFOptimizer,P0,args=(RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest),xtol=1e-100,ftol=1e-100,maxiter=1e100)

                return P.reshape(shape(P0))


        def SDPHDFLocate(self, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
               
                c       = 3e08                              
                if RN_RSS==None:
                        
                        RN_ToA = cvxm.matrix(RN_ToA)
                        ToA = cvxm.matrix(ToA)
                        RoA     = c*ToA
                        RoAStd  = c*ToAStd
                        RoAStd = cvxm.matrix(RoAStd)
                        mtoa,ntoa=cvxm.size(RN_ToA)

                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        RN_TDoA = cvxm.matrix(RN_TDoA)
                        RN_TDoA2 = cvxm.matrix(RN_TDoA2)
                        TDoA = cvxm.matrix(TDoA)
                        RDoA     = c*TDoA
                        RDoAStd=cvxm.matrix(c*TDoAStd)
                        mtdoa,ntdoa=cvxm.size(RN_TDoA)
                
                        Im = cvxm.eye(mtoa)
                        Y=cvxm.optvar('Y',mtoa+1,mtoa+1)
                        t=cvxm.optvar('t',ntoa+ntdoa,1)
                        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
                        prob.constr.append(Y>=0)
                        prob.constr.append(Y[mtoa,mtoa]==1)
                        for i in range(ntoa):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_ToA[:,i])],[-RN_ToA[:,i], cvxm.transpose(RN_ToA[:,i])*RN_ToA[:,i]]])
                            prob.constr.append(-t[i]<(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                            prob.constr.append(t[i]>(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                        for i in range(ntdoa):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_TDoA[:,i])],[-RN_TDoA[:,i], cvxm.transpose(RN_TDoA[:,i])*RN_TDoA[:,i]]])
                            X1=cvxm.matrix([[Im, -cvxm.transpose(RN_TDoA2[:,i])],[-RN_TDoA2[:,i], cvxm.transpose(RN_TDoA2[:,i])*RN_TDoA2[:,i]]])
                            prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
                            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
                        prob.solve()
                        Pval=Y.value
                        X_cvx=Pval[:2,-1]

                        return X_cvx       

                elif RN_ToA==None:

                        RN_RSS=cvxm.matrix(RN_RSS)
                        RSS=cvxm.matrix(RSS.T)
                        RSSnp=cvxm.matrix(RSSnp.T)
                        RSSStd=cvxm.matrix(RSSStd.T)
                        PL0=cvxm.matrix(PL0)
                        mrss,nrss=cvxm.size(RN_RSS)
                
                        Si = array([(1/d0**2)*10**((RSS[0,0]-PL0[0,0])/(5.0*RSSnp[0,0])),(1/d0**2)*10**((RSS[0,1]-PL0[1,0])/(5.0*RSSnp[0,1])),(1/d0**2)*10**((RSS[0,2]-PL0[2,0])/(5.0*RSSnp[0,2])),(1/d0**2)*10**((RSS[0,3]-PL0[3,0])/(5.0*RSSnp[0,3]))])
                
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        RN_TDoA = cvxm.matrix(RN_TDoA)
                        RN_TDoA2 = cvxm.matrix(RN_TDoA2)
                        TDoA = cvxm.matrix(TDoA)
                        RDoA     = c*TDoA
                        RDoAStd=cvxm.matrix(c*TDoAStd)
                        mtdoa,ntdoa=cvxm.size(RN_TDoA)
                
                        Im = cvxm.eye(mrss)
                        Y=cvxm.optvar('Y',mrss+1,mrss+1)
                        t=cvxm.optvar('t',nrss+ntdoa,1)
                        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
                        prob.constr.append(Y>=0)
                        prob.constr.append(Y[mrss,mrss]==1)
                        for i in range(nrss):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_RSS[:,i])],[-RN_RSS[:,i], cvxm.transpose(RN_RSS[:,i])*RN_RSS[:,i]]])
                            prob.constr.append(-RSSStd[0,i]*t[i]<Si[i]*cvxm.trace(X0*Y)-1)
                            prob.constr.append(RSSStd[0,i]*t[i]>Si[i]*cvxm.trace(X0*Y)-1)
                        for i in range(ntdoa):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_TDoA[:,i])],[-RN_TDoA[:,i], cvxm.transpose(RN_TDoA[:,i])*RN_TDoA[:,i]]])
                            X1=cvxm.matrix([[Im, -cvxm.transpose(RN_TDoA2[:,i])],[-RN_TDoA2[:,i], cvxm.transpose(RN_TDoA2[:,i])*RN_TDoA2[:,i]]])
                            prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
                            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
                        prob.solve()
                        Pval=Y.value
                        X_cvx=Pval[:2,-1]

                        return X_cvx
                                        
                elif RN_TDoA==None:
                        
                        RN_RSS=cvxm.matrix(RN_RSS)
                        RSS=cvxm.matrix(RSS.T)
                        RSSnp=cvxm.matrix(RSSnp.T)
                        RSSStd=cvxm.matrix(RSSStd.T)
                        PL0=cvxm.matrix(PL0)
                        mrss,nrss=cvxm.size(RN_RSS)
                
                        Si = array([(1/d0**2)*10**((RSS[0,0]-PL0[0,0])/(5.0*RSSnp[0,0])),(1/d0**2)*10**((RSS[0,1]-PL0[1,0])/(5.0*RSSnp[0,1])),(1/d0**2)*10**((RSS[0,2]-PL0[2,0])/(5.0*RSSnp[0,2])),(1/d0**2)*10**((RSS[0,3]-PL0[3,0])/(5.0*RSSnp[0,3]))])

                        RN_ToA = cvxm.matrix(RN_ToA)
                        ToA = cvxm.matrix(ToA)
                        RoA     = c*ToA
                        RoAStd  = c*ToAStd
                        RoAStd = cvxm.matrix(RoAStd)
                        mtoa,ntoa=cvxm.size(RN_ToA)

                        Im = cvxm.eye(mrss)
                        Y=cvxm.optvar('Y',mrss+1,mrss+1)
                        t=cvxm.optvar('t',nrss+ntoa,1)
                        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
                        prob.constr.append(Y>=0)
                        prob.constr.append(Y[mrss,mrss]==1)

                        for i in range(nrss):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_RSS[:,i])],[-RN_RSS[:,i], cvxm.transpose(RN_RSS[:,i])*RN_RSS[:,i]]])
                            prob.constr.append(-RSSStd[0,i]*t[i]<Si[i]*cvxm.trace(X0*Y)-1)
                            prob.constr.append(RSSStd[0,i]*t[i]>Si[i]*cvxm.trace(X0*Y)-1)
                        for i in range(ntoa):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_ToA[:,i])],[-RN_ToA[:,i], cvxm.transpose(RN_ToA[:,i])*RN_ToA[:,i]]])
                            prob.constr.append(-t[i]<(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                            prob.constr.append(t[i]>(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                        
                        prob.solve()
                        Pval=Y.value
                        X_cvx=Pval[:2,-1]

                        return X_cvx

                else:
                        RN_RSS=cvxm.matrix(RN_RSS)
                        RSS=cvxm.matrix(RSS.T)
                        RSSnp=cvxm.matrix(RSSnp.T)
                        RSSStd=cvxm.matrix(RSSStd.T)
                        PL0=cvxm.matrix(PL0)
                        mrss,nrss=cvxm.size(RN_RSS)
                
                        Si = array([(1/d0**2)*10**((RSS[0,0]-PL0[0,0])/(5.0*RSSnp[0,0])),(1/d0**2)*10**((RSS[0,1]-PL0[1,0])/(5.0*RSSnp[0,1])),(1/d0**2)*10**((RSS[0,2]-PL0[2,0])/(5.0*RSSnp[0,2])),(1/d0**2)*10**((RSS[0,3]-PL0[3,0])/(5.0*RSSnp[0,3]))])

                        RN_ToA = cvxm.matrix(RN_ToA)
                        ToA = cvxm.matrix(ToA)
                        RoA     = c*ToA
                        RoAStd  = c*ToAStd
                        RoAStd = cvxm.matrix(RoAStd)
                        mtoa,ntoa=cvxm.size(RN_ToA)

                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        RN_TDoA = cvxm.matrix(RN_TDoA)
                        RN_TDoA2 = cvxm.matrix(RN_TDoA2)
                        TDoA = cvxm.matrix(TDoA)
                        RDoA     = c*TDoA
                        RDoAStd=cvxm.matrix(c*TDoAStd)
                        mtdoa,ntdoa=cvxm.size(RN_TDoA)

                        Im = cvxm.eye(mrss)
                        Y=cvxm.optvar('Y',mrss+1,mrss+1)
                        t=cvxm.optvar('t',nrss+ntoa+ntdoa,1)
                        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
                        prob.constr.append(Y>=0)
                        prob.constr.append(Y[mrss,mrss]==1)
                        
                        for i in range(nrss):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_RSS[:,i])],[-RN_RSS[:,i], cvxm.transpose(RN_RSS[:,i])*RN_RSS[:,i]]])
                            prob.constr.append(-RSSStd[0,i]*t[i]<Si[i]*cvxm.trace(X0*Y)-1)
                            prob.constr.append(RSSStd[0,i]*t[i]>Si[i]*cvxm.trace(X0*Y)-1)
                        for i in range(ntoa):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_ToA[:,i])],[-RN_ToA[:,i], cvxm.transpose(RN_ToA[:,i])*RN_ToA[:,i]]])
                            prob.constr.append(-t[i]<(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                            prob.constr.append(t[i]>(cvxm.trace(X0*Y)-RoA[i]**2)*(1/RoAStd[i]))
                        for i in range(ntdoa):
                            X0=cvxm.matrix([[Im, -cvxm.transpose(RN_TDoA[:,i])],[-RN_TDoA[:,i], cvxm.transpose(RN_TDoA[:,i])*RN_TDoA[:,i]]])
                            X1=cvxm.matrix([[Im, -cvxm.transpose(RN_TDoA2[:,i])],[-RN_TDoA2[:,i], cvxm.transpose(RN_TDoA2[:,i])*RN_TDoA2[:,i]]])
                            prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
                            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
                        
                        prob.solve()
                        Pval=Y.value
                        X_cvx=Pval[:2,-1]

                        return X_cvx

        def CRBHDFLocate(self, P, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, ToAStd, TDoAStd, PL0, d0, RSSnp, RSSStd, Rest):
                """
                This applies LS approximation to get position P.
                Return P
                """
               
                CRBL=CRBLocation(None)
               
                if RN_RSS==None:
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]
                        #RN_TDoA2= RN_ToA[:,0:1]*ones((1,RNnum_TDoA))
                        return CRBL.CRB_TOA_TDOA_fim(P, RN_ToA, RN_TDoA, RN_TDoA2,ToAStd, TDoAStd)       

                elif RN_ToA==None:
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_RSS[:,0:1]*ones((1,RNnum_TDoA))
                        return CRBL.CRB_RSS_TDOA_fim(P, RN_RSS, RN_TDoA, RN_TDoA2, RSSnp, RSSStd, TDoAStd)
                                        
                elif RN_TDoA==None:
                        return CRBL.CRB_RSS_TOA_fim(P, RN_RSS, RN_ToA, RSSnp, RSSStd, ToAStd)

                else:
                        shRN_TDoA    = shape(RN_TDoA)                                    
                        RNnum_TDoA   = shRN_TDoA[1]                                                                             
                        #RN_TDoA2= RN_TDoA[:,0:1]*ones((1,RNnum_TDoA))
                        return CRBL.CRB_RSS_TOA_TDOA_fim(P, RN_RSS, RN_ToA, RN_TDoA, RN_TDoA2, RSSnp, RSSStd, ToAStd, TDoAStd)
                
                        
