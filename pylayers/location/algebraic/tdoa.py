import os
import numpy as np
import scipy as sp
from scipy import optimize
import numpy.linalg as la
#import cvxmod as cvxm
#import cvxopt as cvxo
from crlb import *

class TDoALocation(object):
    """
    A TDoALocation contains:
    1- 2 sets of RadioNodes (RN1 and RN2) with associated position accuracies (RN1QoS and RN2QoS),
    2- a set of TDoAs measurements (TDoA) with associated STD (TDoAStd)

    This class manages the TDoA-based localization techniques.

    MEMBERS:

        RN1  : an arrray that defines the set of first side Radio nodes implied in localization (coordiantes in meters)
            : shape(RN1)= (2 or 3,RNnum)
        RN1QoS  : an array that defines the precision of positions of RN1 (std in meters)
            : shape(RN1QoS)= (2 or 3, RNnum)
        RN2     : an array that defines the set of second side Radio nodes implied in localization (coordiantes in meters)
            : shape(RN2)= (2 or 3,RNnum)
        RN2QoS  : An Array that defines the precision of positions of RN2 (std in meters)
            : shape(RN2QoS)= (2 or 3, RNnum)
        TDoA    : A measurement vector of TDoA associated to RN (TDoA values in seconds)
            : shape(TDoA)= (RNnum,1)
        TDoAStd : Associated STD of TDoA (std in seconds)
            : shape(TDoAStd)= (RNnum,1)



    Provided Methods:

        info()                      : Display information about the scenario

        LSTDoALocate(RN, TDoA)              : Applies Least Square approximation and estimate position
        WLSTDoALocate(RN, TDoA, TDoAStd)        : Applies Weighted Least Square approximation and estimate position
        TSTDoALocation(P0, RN, TDoA, TDoAStd, Niter)    : Applies Taylor Series method and estimate position after Niter iterations
        TDoAOptimizer(RN, TDoA, TDoAStd)        : Defines the function to be optimized
        MLTDoALocate(P0, RN, TDoA, TDoAStd)         : Optimize TDoAOptimizer() and estimate Position (P0:initial guess)

        CRBTDoALocate(self, P, RN, TDoA, TDoAStd)       : Compute the CRB in P for the given scenario
    """

    def __init__(self,RN1, RN2, TDoA, TDoAStd):
        self.RN1  = RN1
        self.RN2  = RN2
        self.TDoA  = TDoA
        self.TDoAStd  = TDoAStd

    def __init__(self,RN1):
        self.RN1    = RN1


    def info(self):
        """
        Dispaly scenario information
        """
        print("First Reference Radio Nodes:\n" + str(self.RN1))
        print("Second Reference Radio Nodes:\n" + str(self.RN2))
        print("Measured TDoA:\n" + str(self.TDoA))
        print("STD of Measured TDoA:\n" + str(self.TDoAStd))

    def LSTDoALocate(self,RN1, RN2, TDoA):
        """
        This applies LS approximation on TDoA to get position P.
        Return P
        """
        shRN    = shape(RN1)                    # shape of RN
        RNnum   = shRN[1]                       # Number of reference nodes
        c       = 3e08                      # Speed of light
        # Construct the vector K (see theory)
        k1      = (np.sum((RN1-RN2)*(RN1-RN2),axis=0)).reshape(RNnum,1)    # first half of K

        RDoA    = c*TDoA                    # Range of arrival (meters)
        RDoA2   = (RDoA*RDoA).reshape(RNnum,1)
        k2      = RDoA2                     # second half of K

        K       = k1-k2

        # Construct the matrix A (see theory)
        A  = np.hstack((RN1.T - RN2.T,RDoA))

        # Apply LS operator
        Pr = 0.5*np.dot(la.inv(np.dot(A.T,A)),np.dot(A.T,K))
        P = Pr[:shRN[0],:]
        # Return the estimated position
        return P

    def TLSTDoALocate(self,RN1, RN2, TDoA, TDoAStd):
        """
        This applies LS approximation on TDoA to get position P.
        Return P
        """
        shRN    = np.shape(RN1)                    # shape of RN
        RNnum   = shRN[1]                       # Number of reference nodes
        c       = 3e08                      # Speed of light
        # Construct the vector K (see theory)
        k1      = (np.sum((RN1-RN2)*(RN1-RN2),axis=0)).reshape(RNnum,1)    # first half of K

        RDoA    = c*TDoA                    # Range of arrival (meters)
        RDoA2   = (RDoA*RDoA).reshape(RNnum,1)
        k2      = RDoA2                     # second half of K

        K       = k1-k2

        # Construct the matrix A (see theory)
        A  = np.hstack((RN1.T - RN2.T,RDoA))
        A2 = np.dot(transpose(A),A)

        [U,S,V] = la.svd(A2)
        J = 1/S
        rA = la.rank(A)
        m,n = np.shape(A)
        f=0
        if np.log10(la.cond(A2))>=c*max(TDoAStd):
            f=f+1
            for i in range(n-rA):
                u = np.where(J==max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T,np.diag(J)),U.T)
        # Apply LS operator
        Pr      = 0.5*np.dot(A2i,np.dot(A.T,K))
        P       = Pr[:shRN[0],:]
        # Return the estimated position
        return P

    def WLSTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        This applies WLS approximation on TDoA assuming TDoAStd to get position P.
        Return P
        """
        shRN    = shape(RN1)                    # shape of RN
        RNnum   = shRN[1]                       # Number of reference nodes
        c       = 3e08
        RDoAStd = c*TDoAStd                         # Speed of light
        # Construct the vector K (see theory)
        k1      = (np.sum((RN1-RN2)*(RN1-RN2),axis=0)).reshape(RNnum,1)    # first half of K

        RDoA    = c*TDoA                    # Range of arrival (meters)
        RDoA2   = (RDoA*RDoA).reshape(RNnum,1)
        k2      = RDoA2                     # second half of K

        K       = k1-k2

        # Construct the matrix A (see theory)
        A = np.hstack((RN1.T - RN2.T,RDoA))

        # Construct the Covariance Matrix
        C = np.diag(RDoAStd[:,0]**2)

        # Apply LS operator
        Pr = 0.5*dot(linalg.inv(dot(A.T,dot(linalg.inv(C),A))),dot(dot(A.T,linalg.inv(C)),K))
        P = Pr[:shRN[0],:]
        # Return the estimated position
        return P

    def TWLSTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        This applies WLS approximation on TDoA assuming TDoAStd to get position P.
        Return P
        """
        shRN    = np.shape(RN1)                    # shape of RN
        RNnum   = shRN[1]                       # Number of reference nodes
        c       = 3e08
        RDoAStd = c*TDoAStd                         # Speed of light
        # Construct the vector K (see theory)
        k1      = (np.sum((RN1-RN2)*(RN1-RN2),axis=0)).reshape(RNnum,1)    # first half of K

        RDoA    = c*TDoA                    # Range of arrival (meters)
        RDoA2   = (RDoA*RDoA).reshape(RNnum,1)
        k2      = RDoA2                     # second half of K

        K       = k1-k2

        # Construct the matrix A (see theory)
        A       = np.hstack((RN1.T - RN2.T,RDoA))

        # Construct the Covariance Matrix
        C       = np.diag(RDoAStd[:,0]**2)     

        A2   = np.dot(A.T,np.dot(linalg.inv(C),A))

        [U,S,V]= la.svd(A2)
        J = 1/S
        rA = np.rank(A)
        m,n = np.shape(A)
        f=0
        if np.log10(cond(A2))>c*np.max(TDoAStd):
            f=f+1
            for i in range(n-rA):
                u = np.where(J==max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T,np.diag(J)),U.T)
        # Apply LS operator
        Pr      = 0.5*np.dot(A2i,np.dot(np.dot(A.T,la.inv(C)),K))
        P       = Pr[:shRN[0],:]
        # Return the estimated position
        return P

    def TSTDoALocation(self, P0, RN1, RN2, TDoA, TDoAStd, Niter):
        '''
        Applies Taylor Series method and estimate position after Niter iterations
        '''

        P       = P0                        # Initialisation of P as equal to intial guess P0
        shRN    = np.shape(RN1)                    # shape of RN
        RNnum   = shRN[1]                       # Number of reference nodes
        c       = 3e08                      # Speed of light
        RDoA    = c*TDoA
        RDoAStd = c*TDoAStd

        for i in arange(Niter):

            # Construct the matrix A (see theory)
            A  = ((np.outer(P,np.ones(RNnum))-
                RN1)/np.sqrt(sum((np.outer(P,np.ones(RNnum))-
                RN1)**2,axis=0))).T-((np.outer(P,np.ones(RNnum))-
                RN2)/np.sqrt(sum((outer(P,ones(RNnum))- RN2)**2,axis=0))).T

            # Construct the Covariance Matrix
            C  = np.diag((RDoAStd[:,0])**2)

            # Construct the vector D (see theory)
            D = RDoA-(np.sqrt((np.sum((np.outer(P,np.ones(RNnum))-
                RN1)**2,axis=0)).reshape(shape(RDoA)))-sqrt((sum((np.outer(P,np.ones(RNnum))- RN2)**2,axis=0)).reshape(shape(RDoA))))

            # construct the vector Delta (see theory)
            Delta = np.dot(la.inv(np.dot(A.T,np.dot(la.inv(C),A))),np.dot(np.dot(A.T,la.inv(C)),D))
            # update P
            P       = P+Delta
        # Return the estimated position
        return P


    def TDoAOptimizer(self, P, RN1, RN2, TDoA, TDoAStd):
        """
        This defines the ML function to be minimized
        """

        shRN = np.shape(RN1)                    # shape of RN
        RNnum = shRN[1]                       # Number of reference nodes
        c = 3e08                      # Speed of light
        RDoA = c*TDoA
        RDoAStd = c*TDoAStd

        # construct the ML function to be minimized 
        RN1mP = RN1 - np.outer(P, np.ones(RNnum))
        mRN1mP = (np.sqrt(np.diag(np.dot(RN1mP.T,RN1mP)))).reshape(RNnum,1)
        RN2mP  = RN2 - np.outer(P,np.ones(RNnum))
        mRN2mP = (np.sqrt(np.diag(np.dot(RN2mP.T,RN2mP)))).reshape(RNnum,1)
        rRDoA = mRN1mP-mRN2mP

        tk = (RDoA-rRDoA)**2/(2*RDoAStd**2)
        uk = tk[:,0]# *(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])
        suk  = uk.sum(axis=0)
        #msuk    = sqrt(dot(suk,suk.T))

        return(suk)
    

    def MLTDoALocate(self, P, P0, RN1, RN2, TDoA, TDoAStd):
        """ 
        Optimization Routine
        """
        P = optimize.fmin(self.TDoAOptimizer,P0,args=(RN1, RN2,TDoA,TDoAStd),xtol=1e-10,ftol=1e-10)
        return P.reshape(np.shape(P0))

    '''def SDPTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        Apply SDP approximation and localization
        """
        RN1 = cvxm.matrix(RN1)
        RN2 = cvxm.matrix(RN2)
        TDoA = cvxm.matrix(TDoA)
        c       = 3e08                      
        RDoA     = c*TDoA
        RDoAStd=cvxm.matrix(c*TDoAStd)
        mtdoa,ntdoa=cvxm.size(RN1)
        Im = cvxm.eye(mtdoa)
        Y=cvxm.optvar('Y',mtdoa+1,mtdoa+1)
        t=cvxm.optvar('t',ntdoa,1)
        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
        prob.constr.append(Y>=0)
        prob.constr.append(Y[mtdoa,mtdoa]==1)
        for i in range(ntdoa):
            X0=cvxm.matrix([[Im, -cvxm.transpose(RN1[:,i])],[-RN1[:,i], cvxm.transpose(RN1[:,i])*RN1[:,i]]])
            X1=cvxm.matrix([[Im, -cvxm.transpose(RN2[:,i])],[-RN2[:,i], cvxm.transpose(RN2[:,i])*RN2[:,i]]])
            prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace(X0*Y)+cvxm.trace(X1*Y)-RDoA[i,0]**2)
            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace(X0*Y)+cvxm.trace(X1*Y)-RDoA[i,0]**2)
        prob.solve()
        Pval=Y.value
        X_cvx=Pval[:2,-1]

        return X_cvx'''

#    def SDPTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
#        """
#        Apply SDP approximation and localization
#        """
#        RN1 = cvxm.matrix(RN1)
#        RN2 = cvxm.matrix(RN2)
#        TDoA = cvxm.matrix(TDoA)
#        c       = 3e08                      
#        RDoA     = c*TDoA
#        RDoAStd=cvxm.matrix(c*TDoAStd)
#        mtdoa,ntdoa=cvxm.size(RN1)
#        Im = cvxm.eye(mtdoa)
#        Y=cvxm.optvar('Y',mtdoa+1,mtdoa+1)
#        t=cvxm.optvar('t',ntdoa,1)
#        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
#        prob.constr.append(Y>=0)
#        prob.constr.append(Y[mtdoa,mtdoa]==1)
#        for i in range(ntdoa):
#            X0=cvxm.matrix([[Im, -cvxm.transpose(RN1[:,i])],[-RN1[:,i], cvxm.transpose(RN1[:,i])*RN1[:,i]]])
#            X1=cvxm.matrix([[Im, -cvxm.transpose(RN2[:,i])],[-RN2[:,i], cvxm.transpose(RN2[:,i])*RN2[:,i]]])
#            '''prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)
#            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace(X0*Y)-cvxm.trace(X1*Y)-RDoA[i,0]**2)'''
#            prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace((X1-X0)*Y)-RDoA[i,0]**2)
#            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace((X1-X0)*Y)-RDoA[i,0]**2)
#        prob.solve()
#        Pval=Y.value
#        X_cvx=Pval[:2,-1]
#
#        return X_cvx


#    def SDPTDoALocate1(self, RN1, RN2, TDoA, TDoAStd):
#        """
#        Apply SDP approximation and localization
#        """
#        RN1 = cvxm.matrix(RN1)
#        RN2 = cvxm.matrix(RN2)
#        TDoA = cvxm.matrix(TDoA)
#        c       = 3e08                      
#        RDoA     = c*TDoA
#        RDoAStd=cvxm.matrix(c*TDoAStd)
#        mtdoa,ntdoa=cvxm.size(RN1)
#        Im = cvxm.eye(mtdoa)
#        Y=cvxm.optvar('Y',mtdoa+1,mtdoa+1)
#        t=cvxm.optvar('t',ntdoa,1)
#        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
#        prob.constr.append(Y>=0)
#        prob.constr.append(Y[mtdoa,mtdoa]==1)
#        for i in range(ntdoa):
#            
#            X0=cvxm.matrix([[Im, -cvxm.transpose(RN1[:,i])],[-RN1[:,i], cvxm.transpose(RN1[:,i])*RN1[:,i]]])
#            X1=cvxm.matrix([[Im, -cvxm.transpose(RN2[:,i])],[-RN2[:,i], cvxm.transpose(RN2[:,i])*RN2[:,i]]])
#            prob.constr.append(-RDoAStd[i,0]*t[i]<cvxm.trace((X1-X0)*Y)-RDoA[i,0]**2)
#            prob.constr.append(RDoAStd[i,0]*t[i]>cvxm.trace((X1-X0)*Y)-RDoA[i,0]**2)
#        prob.solve()
#        Pval=Y.value
#        X_cvx=Pval[:2,-1]
#
#        return X_cvx

    def CRBTDoALocate(self, P, RN1, RN2, TDoAStd):
        """
        Compute the CRB in P for the given scenario
        """
        '''c       = 3e08
        shP     = shape(P)
        shRN    = shape(RN1)
        RNnum   = shRN[1]
        RDoAStd = c*TDoAStd

        RN1mP   = outer(P,ones(RNnum))- RN1
        mRN1mP  = (sqrt(diag(dot(RN1mP.T,RN1mP)))).reshape(RNnum,1)
        RN2mP   = outer(P,ones(RNnum))- RN2
        mRN2mP  = (sqrt(diag(dot(RN2mP.T,RN2mP)))).reshape(RNnum,1)

        num     = sum(2/(RDoAStd[:,0]**2)*(1-sum((RN1mP/mRN1mP[:,0])*(RN2mP/mRN2mP[:,0]),axis=0)),axis=0)       # the numerator of the CRLB


        div1    = sum((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])**2,axis=1).reshape(shP)
        don1    = div1.prod(axis=0)[0]      # first term of the doniminator

        div2    = sum(prod((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0]),axis=0),axis=0)
        don2    = div2**2               # second term of the doniminator

        CRB     = num/(don1-don2)           # the CRB'''
        crlb=CRBLocation(RN1)
        #CRB=crlb.CRB_TDOA_fim(P, RN1, RN2, TDoAStd)
        CRB=crlb.Angle_TDOA(P, RN1, RN2, TDoAStd)

        return sqrt(abs(CRB[0]))



