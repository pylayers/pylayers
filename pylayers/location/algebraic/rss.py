import os
import numpy as np
import scipy as sp
from scipy import optimize
import numpy.linalg as la
#import cvxmod as cvxm
#import cvxopt as cvxo
from string import *
from crlb import *

class RSSLocation(object):
    """
    A RSSALocation contains:
      1- a set of RadioNodes (RN) with associated position accuracies (RNQoS),
      2- a set of RSSs measurements (RSS) with associated STD of shadowing(RSSStd) and associated propagation constant (RSSnp) 

    This class manages the RSS-based localization techniques.
    MEMBERS:

    RN      : An Array that defines the Radio nodes implied in localization (coordiantes in meters)
                        : shape(RN)= (2 or 3,RNnum)
    RNQoS   : An Array that defines the precision of positions of RN (std in meters)
                        : shape(RNQoS)= (2 or 3, RNnum)
    RSS     : A measurement vector of RSS associated to RN (RSS values in dB)
                        : shape(RSS)= (RNnum,1)
    RSSStd  : Associated STD of shadowing (std in dB)
                        : shape(RSSStd)= (RNnum,1)
    RSSnp   : Associated propagation constant
                        : shape(RSSnp)= (RNnum,1)
    lamda   : Associated wavelength (meters)


    Provided Methods:
     info() : Display information about the scenario

     getPL0(lamda, d0)                           : Compute PL0
     getPLmean(RN, PL0, d0, RSSnp)               : Compute PL mean
     getPL(RN, lamda, d0, RSSnp, RSSStd)         : Compute PL

     getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    : Compute Ranges using "Rest" estimator from RSS
     getRangeStd(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)  : Compute Ranges std associated to "Rest" estimator
     LSRSSLocate(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    : Applies Least Square approximation and estimate position
     WLSRSSLocate(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)   : Applies Weighted Least Square approximation and estimate position
     IRSSOptimizer(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)  : Defines the function to be optimized (indirect estimator)
     MLIRSSLocate(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)   : Optimize IRSSOptimizer() and estimate Position (P0:initial guess)
     DRSSOptimizer(RN, PL0, d0, RSS, RSSnp, RSSStd)        : Defines the function to be optimized (direct estimator)
     MLDRSSLocate(RN, PL0, d0, RSS, RSSnp, RSSStd)         : Optimize DRSSOptimizer() and estimate Position (P0:initial guess)

     CRBRSSLocate(self, P, RN, PL0, d0, RSS, RSSnp, RSSStd): This computes CRB of RSS positioning
        """

    def __init__(self, RN, PL0, d0, RSS, RSSnp, RSSStd):
        self.RN         = RN
        self.PL0        = PL0
        self.d0         = d0
        self.RSS        = RSS
        self.RSSnp      = RSSnp
        self.RSSStd     = RSSStd

    def __init__(self, RN, PL0, d0, RSSnp, RSSStd):
        self.RN         = RN
        self.PL0        = PL0
        self.d0         = d0
        self.RSSnp      = RSSnp
        self.RSSStd     = RSSStd
    def __init__(self, RN):
        self.RN         = RN

    def info(self):
        """
        Display scenario information
        """
        print "Reference Radio Nodes:\n", self.RN
        print "References distances:\n", self.d0
        print "RSSI at d0:\n", self.lamda
        print "Measured RSS:\n", self.RSS
        print "Propagation constants:\n", self.RSSnp
        print "STD of Measured RSS shadowing:\n", self.RSSStd

    def getPL0(self, lamda, d0):
       """ Compute PL0

       Parameters
       ----------

       """
       return  20*np.log10(4*np.pi*d0/lamda)

    def getPLmean(self, RN, P, PL0, d0, RSSnp):
       """ Compute PL mean
       """
       shRN            = shape(RN)
       RNmP            = (np.sqrt(sum((RN-P)**2,axis=0))).reshape((shRN[1],1))                    # distance between AN and P
       return          PL0-10*RSSnp*np.log10(RNmP/d0)

    def getPL(self, RN, P, PL0, d0, RSSnp, RSSStd):
        """ Compute PL

        Parameters
        ----------

        """
        PLmean      = self.getPLmean(RN, P, PL0, d0, RSSnp)
        shPLmean    = np.shape(PLmean)
        Xrand       = RSSStd*sp.randn(shPLmean[0],shPLmean[1])

        return      PLmean+Xrand

    def getRange(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """ Compute Ranges using "Rest" estimator from RSS

        Parameters
        ----------

        """
        S           = -(np.log(10)/10)* RSSStd/RSSnp            # STD of ranges distribution
        M           =  (np.log(10)/10)*(PL0-RSS)/RSSnp + np.log(d0)    # Mean of ranges distribution


        if lower(Rest)      == 'mode':
            return      np.exp(M-S**2)
        elif lower(Rest)    == 'median':
            return      np.exp(M)
        elif lower(Rest)    == 'mean':
            return      np.exp(M+0.5*S**2)
        else:
            return      np.exp(M)
            print       "No \"%s\" defined estimator" %Rest

    def getRangeStd(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """
        Compute Ranges std associated to "Rest" estimator
        """
        S  = -(np.log(10)/10)* RSSStd/RSSnp                    # STD of ranges distribution
        M  = (np.log(10)/10)*(PL0-RSS)/RSSnp + np.log(d0)         # Mean of ranges distribution
        if lower(Rest)      == 'mode':
            return      np.sqrt((np.exp(2*M-2*S**2))*(-np.exp(-S**2)+1))
        elif lower(Rest)    == 'median':
            return      np.sqrt((np.exp(2*M+S**2))*(np.exp(S**2)-1))
        elif lower(Rest)    == 'mean':
            return      np.sqrt((np.exp(2*M+3*S**2))*(np.exp(S**2)-1))
        else:
            return      np.sqrt((np.exp(2*M+S**2))*(np.exp(S**2)-1))
            print       "No \"%s\" defined estimator" %Rest

    def LSRSSLocate(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """
        This applies LS approximation on RSS based ranges to get position P.
        Return P
        """
        shRN    = shape(RN)                         # shape of RN
        RNnum   = shRN[1]                           # Number of reference nodes
        c       = 3e08                          # Speed of light
        # Construct the vector K (see theory)
        RN2     = (np.sum(RN*RN,axis=0)).reshape(RNnum,1)
        k1      = RN2[1:RNnum,:]-RN2[0,0]                   # first half of K
        RoA     = self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    # RSS based Ranges (meters)
        RoA2    = (RoA*RoA).reshape(RNnum,1)
        k2      = RoA2[0,0]-RoA2[1:RNnum,:]                 # second half of K
        K       = k1+k2
        # Construct the matrix A (see theory)
        A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])

        # Apply LS operator
        P       = 0.5*np.dot(la.inv(np.dot(A.T,A)),np.dot(A.T,K))
        # Return the estimated position
        return P

    def TLSRSSLocate(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """
        This applies TLS approximation on RSS based ranges to get position P.
        Return P
        """
        shRN    = np.shape(RN)                         # shape of RN
        RNnum   = shRN[1]                           # Number of reference nodes
        c       = 3e08                          # Speed of light
        # Construct the vector K (see theory)
        RN2     = (np.sum(RN*RN,axis=0)).reshape(RNnum,1)
        k1      = RN2[1:RNnum,:]-RN2[0,0]                   # first half of K
        RoA     = self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    # RSS based Ranges (meters)
        RoA2    = (RoA*RoA).reshape(RNnum,1)
        k2      = RoA2[0,0]-RoA2[1:RNnum,:]                 # second half of K
        K       = k1+k2
        # Construct the matrix A (see theory)
        A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])

        A2   = np.dot(np.transpose(A),A)

        [U,S,V]=la.svd(A2)
        J = 1/S
        rA =  la.rank(A)
        m,n = np.shape(A)
        f=0
        if np.log10(cond(A2))>=max(self.getRangeStd(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
            f=f+1
            for i in range(n-rA):
                u = where(J==max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T,la.diag(J)),U.T)

        P       = 0.5*np.dot(A2i,np.dot(A.T,K))
        # Return the estimated position
        return P

    def WLSRSSLocate(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """ applies WLS approximation on RSS assuming RSSStd to get position P.

        Returns
        -------

        P : estimated position

        """
        shRN    = np.shape(RN)                         # shape of RN
        RNnum   = shRN[1]                           # Number of reference nodes
        # Construct the vector K (see theory)
        RN2     = (np.sum(RN*RN,axis=0)).reshape(RNnum,1)
        k1      = RN2[1:RNnum,:]-RN2[0,0]                   # first half of K
        RoA     = self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    # RSS based Ranges (meters)
        RoAStd  = self.getRangeStd(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        RoA2    = (RoA*RoA).reshape(RNnum,1)
        k2      = RoA2[0,0]-RoA2[1:RNnum,:]             # second half of K
        K       = k1+k2
        # Construct the matrix A (see theory)
        A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])
        # Construct the Covariance Matrix
        C       = la.diag((RoAStd[1:RNnum,0])**2)
        # Apply LS operator
        P       = 0.5*np.dot(la.inv(np.dot(A.T,np.dot(la.inv(C),A))),
                             np.dot(np.dot(A.T,la.inv(C)),K))

        # Return the estimated position
        return P

    def TWLSRSSLocate(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """
        This applies WLS approximation on RSS assuming RSSStd to get position P.
        Return P
        """
        shRN    = np.shape(RN)                         # shape of RN
        RNnum   = shRN[1]                           # Number of reference nodes
        # Construct the vector K (see theory)
        RN2     = (np.sum(RN*RN,axis=0)).reshape(RNnum,1)
        k1      = RN2[1:RNnum,:]-RN2[0,0]                   # first half of K
        RoA     = self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    # RSS based Ranges (meters)
        RoAStd  = self.getRangeStd(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        RoA2    = (RoA*RoA).reshape(RNnum,1)
        k2      = RoA2[0,0]-RoA2[1:RNnum,:]             # second half of K
        K       = k1+k2
        # Construct the matrix A (see theory)
        A       = RN[:,1:RNnum].T - RN[:,0].reshape(1,shRN[0])
        # Construct the Covariance Matrix
        C       = la.diag((RoAStd[1:RNnum,0])**2)
        
        A2   = np.dot(A.T,np.dot(la.inv(C),A))

        [U,S,V] = la.svd(A2)
        J = 1/S
        rA = la.rank(A)
        m,n = np.shape(A)
        f=0
        if np.log10(cond(A2))>=max(self.getRangeStd(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)):
            f=f+1
            for i in range(n-rA):
                u = np.where(J==max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T,la.diag(J)),U.T)
        P   = 0.5*np.dot(A2i,np.dot(np.dot(A.T,la.inv(C)),K))
        return P

    def IRSSOptimizer(self, P, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """
        This defines the ML function to be minimized in the Indirect case
        """

        shRN    = np.shape(RN)                     # shape of RN
        RNnum   = shRN[1]                       # Number of reference nodes
        RoA     = self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)    # RSS based Ranges (meters)
        RoAStd  = self.getRangeStd(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)

        # construct the ML function to be minimized
        RNmP    = RN - np.outer(P,np.ones(RNnum))
        mRNmP   = (np.sqrt(diag(np.dot(RNmP.T,RNmP)))).reshape(RNnum,1)
        tk      = (RoA-mRNmP)**2
        uk      = tk/(2*RoAStd**2)+np.log(np.sqrt(2*np.pi)*RoAStd)
        suk     = uk.sum(axis=0)
        msuk    = suk

        return(msuk)


    def MLIRSSLocate(self, P, P0, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
        """ Optimization Routine of the indirect case
        """

        P = optimize.fmin(self.IRSSOptimizer,P0,args=(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest),xtol=1e-10,ftol=1e-10)

        return P.reshape(shape(P0))

    def DRSSOptimizer(self, P, RN, PL0, d0, RSS, RSSnp, RSSStd):
        """
        This defines the ML function to be minimized in the direct case
        """

        shRN    = shape(RN)             # shape of RN
        RNnum   = shRN[1]               # Number of reference nodes
        S       = -(np.log(10)/10)* RSSStd/RSSnp                    # STD of ranges distribution
        M       = (np.log(10)/10)*(PL0-RSS)/RSSnp + np.log(d0)         # Mean of ranges distribution

        # construct the ML function to be minimized
        RNmP    = RN - np.outer(P,ones(RNnum))
        mRNmP   = (np.sqrt(diag(dot(RNmP.T,RNmP)))).reshape(RNnum,1)
        tk      = (M-S**2-np.log(mRNmP))**2
        uk      = tk/(2*S**2)
        suk     = uk.sum(axis=0)

        
        return(suk)

    def DRSSOptimizer1(self, P, RN, PL0, d0, RSS, RSSnp, RSSStd):
        """
        This defines the ML function to be minimized in the direct case
        """

        shRN    = shape(RN)             # shape of RN
        RNnum   = shRN[1]               # Number of reference nodes
        
        S           = -(np.log(10)/10)* RSSStd/RSSnp                    # STD of ranges distribution
        M           = (np.log(10)/10)*(PL0-RSS)/RSSnp + np.log(d0)         # Mean of ranges distribution
        

        # construct the ML function to be minimized     
        RNmP    = RN - np.outer(P,ones(RNnum))
        mRNmP   = (np.sqrt(diag(dot(RNmP.T,RNmP)))).reshape(RNnum,1)
        tk      = (RSS-PL0-10*RSSnp*np.log10(mRNmP/d0))**2
        uk      = tk/(2*RSSStd**2)
        suk     = uk.sum(axis=0)
        
        
        return(suk)

    def MLDRSSLocate(self, P, P0, RN, PL0, d0, RSS, RSSnp, RSSStd):
        """ 
        Optimization Routine of the direct case
        """

        P       = optimize.fmin(self.DRSSOptimizer,P0,args=(RN, PL0, d0, RSS, RSSnp, RSSStd),xtol=1e-10,ftol=1e-10)

        return P.reshape(shape(P0))

    '''def SDPRSSLocate(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):

        RoA=self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)

        RN=cvxm.matrix(RN)
        RSS=cvxm.matrix(RSS.T)
        RSSnp=cvxm.matrix(RSSnp.T)
        RSSStd=cvxm.matrix(RSSStd.T)
        RoA=cvxm.matrix(RoA.T)
        mrss,nrss=cvxm.size(RN)
        
        Si = array([(RSS[0,0]*10**(RSSStd[0,0]/10.0))/(RoA[0,0]**RSSnp[0,0]),(RSS[0,1]*10**(RSSStd[0,1]/10.0))/(RoA[0,1]**RSSnp[0,1]),(RSS[0,2]*10**(RSSStd[0,2]/10.0))/(RoA[0,2]**RSSnp[0,2]),(RSS[0,3]*10**(RSSStd[0,3]/10.0))/(RoA[0,3]**RSSnp[0,3])])
        qi = cvxm.matrix((Si/RSS)**(2/RSSnp[0,0]))
        Im = cvxm.eye(mrss)
        Y=cvxm.optvar('Y',mrss+1,mrss+1)
        t=cvxm.optvar('t',nrss,1)

        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
        prob.constr.append(Y>=0)
        prob.constr.append(Y[mrss,mrss]==1)
        for i in range(nrss):
        X0=cvxm.matrix([[Im, -cvxm.transpose(RN[:,i])],[-RN[:,i], cvxm.transpose(RN[:,i])*RN[:,i]]])
        prob.constr.append(-t[i]<qi[i]*cvxm.trace(X0*Y)-1)
        prob.constr.append(t[i]>qi[i]*cvxm.trace(X0*Y)-1)
       
        prob.solve()
        Pval=Y.value
        X_cvx=Pval[:2,-1]
        return X_cvx'''

#    def SDPRSSLocate(self, RN, PL0, d0, RSS, RSSnp, RSSStd, Rest):
#
#        RoA=self.getRange(RN, PL0, d0, RSS, RSSnp, RSSStd, Rest)
#
#        RN=cvxm.matrix(RN)
#        RSS=cvxm.matrix(RSS)
#        RSSnp=cvxm.matrix(RSSnp)
#        RSSStd=cvxm.matrix(RSSStd)
#        PL0=cvxm.matrix(PL0)
#        RoA=cvxm.matrix(RoA)
#        mrss,nrss=cvxm.size(RN)
#        Si = array([(1/d0**2)*10**((RSS[0,0]-PL0[0,0])/(5.0*RSSnp[0,0])),(1/d0**2)*10**((RSS[1,0]-PL0[1,0])/(5.0*RSSnp[1,0])),(1/d0**2)*10**((RSS[2,0]-PL0[2,0])/(5.0*RSSnp[2,0])),(1/d0**2)*10**((RSS[3,0]-PL0[3,0])/(5.0*RSSnp[3,0]))])
#        #Si = array([(1/d0**2)*10**(-(RSS[0,0]-PL0[0,0])/(5.0*RSSnp[0,0])),(1/d0**2)*10**(-(RSS[0,1]-PL0[1,0])/(5.0*RSSnp[0,1])),(1/d0**2)*10**(-(RSS[0,2]-PL0[2,0])/(5.0*RSSnp[0,2])),(1/d0**2)*10**(-(RSS[0,3]-PL0[3,0])/(5.0*RSSnp[0,3]))])
#        Im = cvxm.eye(mrss)
#        Y=cvxm.optvar('Y',mrss+1,mrss+1)
#        t=cvxm.optvar('t',nrss,1)
#        prob=cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
#        prob.constr.append(Y>=0)
#        prob.constr.append(Y[mrss,mrss]==1)
#        for i in range(nrss):
#            X0 = cvxm.matrix([[Im, -cvxm.transpose(RN[:,i])],[-RN[:,i], cvxm.transpose(RN[:,i])*RN[:,i]]])
#            prob.constr.append(-RSSStd[i,0]*t[i]<Si[i]*cvxm.trace(X0*Y)-1)
#            prob.constr.append(RSSStd[i,0]*t[i]>Si[i]*cvxm.trace(X0*Y)-1)
#
#        prob.solve()
#        Pval=Y.value
#        X_cvx=Pval[:2,-1]
#        return X_cvx
#

    def CRBRSSLocate(self, P, RN, PL0, d0, RSSnp, RSSStd):
        """
        This computes CRB of RSS positioning
        """
        shP     = shape(P)
        shRN    = shape(RN)
        RNnum   = shRN[1]

        S       = (np.log(10)/10)* RSSStd/RSSnp

        RNmP    = RN - np.outer(P,ones(RNnum))
        mRNmP   = (np.sqrt(diag(dot(RNmP.T,RNmP))))

        num     = sum((1+S**2)/((S**2)*mRNmP**2),axis=0)[0]     # the numerator of the CRLB

        div1    = sum(((1+S[:,0]**2)*RNmP**2)/((S[:,0]**2)*mRNmP**4),axis=1).reshape(shP)
        don1    = div1.prod(axis=0)[0]      # first term of the doniminator

        div2    = sum(((1+S[:,0]**2)*(RNmP.prod(axis=0)))/((S[:,0]**2)*mRNmP**4),axis=0)
        don2    = div2**2               # second term of the doniminator

        CRB     = num/(don1-don2)           # the CRB'''
        crlb=CRBLocation(RN)
        CRB=crlb.CRB_RSS_fim(P, RN, RSSnp, RSSStd)
        return np.sqrt(CRB)
