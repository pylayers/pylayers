import os
import numpy as np
import scipy as sp
from scipy import optimize
import string as stri
import numpy.linalg as la


class CRBLocation(object):
    """
    A CRBLocation contains:
    1- a set of RadioNodes (RN) with associated position accuracies (RNQoS),
    2- a set of measurements (RSS, TOA, TDOA) with associated accuracies.

    This class manages the CRB techniques.
    MEMBERS:

    RN      : array that defines the Radio nodes implied in localization (coordiantes in meters)
                    : shape(RN)= (2 or 3,RNnum)
    RNQoS   : array that defines the precision of positions of RN (std in meters)
                    : shape(RNQoS)= (2 or 3, RNnum)
    param   : a set of parameters depending on the type of measurement [list]. 


    """

    def __init__(self, RN):
        self.RN = RN



    def info(self):
        """
            Display scenario information
        """
        print "Reference Radio Nodes:\n", self.RN
        print "parameters:\n", self.param

    def FIM_RSS(self, P, RN_RSS, RSSnp, RSSStd):
        """  FIM in P of RSS positioning
        """
        shP     = np.shape(P)
        shRN    = np.shape(RN_RSS)
        RNnum   = np.shRN[1]

        S  = (np.log(10)/10)* RSSStd/RSSnp

        RNmP    = RN_RSS - np.outer(P,np.ones(RNnum))
        mRNmP   = (np.sqrt(np.diag(np.dot(RNmP.T,RNmP))))
        j11     = np.sum(((1+S[:,0]**2)*RNmP[0,:]**2)/((S[:,0]**2)*mRNmP**4),axis=0)
        j22     = np.sum(((1+S[:,0]**2)*RNmP[1,:]**2)/((S[:,0]**2)*mRNmP**4),axis=0)
        j12=j21 = np.sum(((1+S[:,0]**2)*(RNmP.prod(axis=0)))/((S[:,0]**2)*mRNmP**4),axis=0)
        FIM  = np.array([[j11,j12],[j21,j22]])
        return FIM

    def FIM_TOA(self, P, RN_TOA, TOAStd):
        """
            Compute the FIM in P for the given scenario
        """
        c       = 3e08
        shP     = np.shape(P)
        shRN    = np.shape(RN_TOA)
        RNnum   = shRN[1]
        RoAStd  = c*TOAStd

        num     = np.sum(1/(RoAStd**2),axis=0)[0]  # the numerator of the CRLB
        RNmP    = RN_TOA - np.outer(P,np.ones(RNnum))
        mRNmP   = (np.sqrt(np.diag(np.dot(RNmP.T,RNmP))))

        j11     = np.sum(RNmP[0,:]**2/((RoAStd[:,0]**2)*mRNmP**2),axis=0)
        j22     = np.sum(RNmP[1,:]**2/((RoAStd[:,0]**2)*mRNmP**2),axis=0)
        j12=j21 = np.sum((RNmP.prod(axis=0))/((RoAStd[:,0]**2)*mRNmP**2),axis=0)
        FIM     = array([[j11,j12],[j21,j22]])

        return FIM

    def FIM_TDOA(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """ Compute the FIM in P for the given scenario

        Parameters
        ----------

        P
        RN1_TDOA
        RN2_TDOA
        TDOAStd

        """
        c  = 0.3
        shP = np.shape(P)
        shRN = np.shape(RN1_TDOA)
        RNnum = shRN[1]
        RDoAStd = c*TDOAStd
        RN1mP = np.outer(P,np.ones(RNnum))- RN1_TDOA
        mRN1mP = (np.sqrt(np.diag(np.dot(RN1mP.T,RN1mP)))).reshape(RNnum,1)
        RN2mP = np.outer(P,np.ones(RNnum))- RN2_TDOA
        mRN2mP = (np.sqrt(np.diag(np.dot(RN2mP.T,RN2mP)))).reshape(RNnum,1)

        j11 = np.sum((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])**2,axis=1)[0]
        j22 = np.sum((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])**2,axis=1)[1]
        j12 = j21 = np.sum(np.prod((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0]),axis=0),axis=0)
        j12a = j21a = np.sum(np.prod((1/RDoAStd[:,0]**2)*(-RN1mP/mRN1mP[:,0]+RN2mP/mRN2mP[:,0]),axis=0),axis=0)
        FIM = np.array([[j11,j12],[j21,j22]])
        return FIM

    def CRB_RSS_fim(self, P, RN_RSS, RSSnp, RSSStd):
        """ computes CRB of RSS positioning as the trace of inv of fim

        Parameters
        ----------

        P
        RN_RSS
        RSSnp
        RSSSstd

        """
        FIM = self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd)

        return np.trace(la.inv(FIM))

    def CRB_TOA_fim(self, P, RN_TOA, TOAStd):
        """ compute the CRB in P for the given scenario  as the trace of inv of fim
        """
        FIM = self.FIM_TOA(P, RN_TOA, TOAStd)

        return np.trace(la.inv(FIM))

    def CRB_TDOA_fim(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """ compute the CRB in P for the given scenario as the trace of inv of fim
        """
        FIM=self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return np.trace(la.inv(FIM))

    def CRB_RSS_TOA_fim(self, P, RN_RSS, RN_TOA, RSSnp, RSSStd, TOAStd):
        """ compute CRB of RSS/TOA positioning as the trace of inv of fim
        """

        FIM=self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd)+self.FIM_TOA(P, RN_TOA, TOAStd)
        return np.trace(la.inv(FIM))

    def CRB_RSS_TDOA_fim(self, P, RN_RSS, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TDOAStd):
        """ computes CRB of RSS/TDOA positioning as the trace of inv of fim
        """
        FIM=self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd)+self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return np.trace(la.inv(FIM))

    def CRB_TOA_TDOA_fim(self, P, RN_TOA, RN1_TDOA, RN2_TDOA,TOAStd, TDOAStd):
        """ computes CRB of TOA/TDOA positioning as the trace of inv of fim
        """
        FIM=self.FIM_TOA(P, RN_TOA, TOAStd)+self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return np.trace(la.inv(FIM))

    def CRB_RSS_TOA_TDOA_fim(self, P, RN_RSS, RN_TOA, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TOAStd, TDOAStd):
        """ computes CRB of RSS/TOA/TDOA positioning as the trace of inv of fim
        """
        FIM=self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd)+self.FIM_TOA(P, RN_TOA, TOAStd)+self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return np.trace(la.inv(FIM))
    def MCRB_RSS_fim(self, L, RN_RSS, RSSnp, RSSStd):
        """
            This computes mean CRB of RSS positioning over the area of length L
        """
        delta = L/50.0
        CRB=[]
        for x in np.arange(0.001,L+0.1,delta):
            for y in np.arange(0.001,L+0.1,delta):
                P=np.array([[x],[y]])
                f1 = self.CRB_RSS_fim(P, RN_RSS, RSSnp, RSSStd)
                if np.isnan(f1)==0 :
                    CRB.append(f1)
        moy = np.mean(CRB)
        return moy

    def MCRB_TOA_fim(self, L, RN_TOA, TOAStd):
        """ computes mean CRB of TOA positioning over the area of length L
        """
        delta = L/50.0
        CRB=[]
        for x in np.arange(0.001,L+0.1,delta): 
            for y in np.arange(0.001,L+0.1,delta): 
                P = np.array([[x],[y]]) 
                f1=self.CRB_TOA_fim(P, RN_TOA, TOAStd)
                if np.isnan(f1)==0 :
                    CRB.append(f1)
        moy = np.mean(CRB)
        return moy

    def MCRB_TDOA_fim(self, L, RN1_TDOA, RN2_TDOA, TDOAStd):
        """ computes mean CRB of TDOA positioning over the area of length L
        """
        delta = L/50.0
        CRB = []
        for x in np.arange(0.001,L+0.1,delta):
            for y in np.arange(0.001,L+0.1,delta):
                P=np.array([[x],[y]])
                f1=self.CRB_TDOA_fim(P, RN1_TDOA, RN2_TDOA, TDOAStd)
                if np.isnan(f1)==0 :
                    CRB.append(f1)
        moy = np.mean(CRB)
        return moy

    def MCRB_RSS_TOA_fim(self, L, RN_RSS, RN_TOA, RSSnp, RSSStd, TOAStd):
        """ computes mean CRB of RSS/TOA positioning over the area of length L
        """
        delta = L/50.0
        CRB=[]
        for x in np.arange(0.001,L+0.1,delta):
            for y in np.arange(0.001,L+0.1,delta):
                P=np.array([[x],[y]])
                f1=self.CRB_RSS_TOA_fim(P, RN_RSS, RN_TOA, RSSnp, RSSStd, TOAStd)
                if not np.isnan(f1):
                    CRB.append(f1)
        moy=np.mean(CRB)
        return moy

    def MCRB_RSS_TDOA_fim(self, L, RN_RSS, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TDOAStd):
        """ This computes mean CRB of RSS/TDOA positioning over the area of length L
        """
        delta = L/50.0
        CRB=[]
        for x in np.arange(0.001,L+0.1,delta): 
            for y in np.arange(0.001,L+0.1,delta):
                P=np.array([[x],[y]])
                f1=self.CRB_RSS_TDOA_fim(P, RN_RSS, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TDOAStd)
                if np.isnan(f1)==0 :
                    CRB.append(f1)

        moy = np.mean(CRB)
        return moy

    def MCRB_TOA_TDOA_fim(self, L, RN_TOA, RN1_TDOA, RN2_TDOA,TOAStd, TDOAStd):
        """ This computes mean CRB of TOA/TDOA positioning over the area of length L
        """
        delta = L/50.0
        CRB=[]
        for x in np.arange(0.001,L+0.1,delta): 
            for y in np.arange(0.001,L+0.1,delta): 
                P = np.array([[x],[y]])
                f1 = self.CRB_TOA_TDOA_fim(P, RN_TOA, RN1_TDOA, RN2_TDOA,TOAStd, TDOAStd)
                if not np.isnan(f1):
                    CRB.append(f1)
    
        moy = np.mean(CRB)
        return moy
    def MCRB_RSS_TOA_TDOA_fim(self, L, RN_RSS, RN_TOA, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TOAStd, TDOAStd):
            """
            This computes mean CRB of RSS/TOA/TDOA positioning over the area of length L
            """

            delta = L/50.0
            CRB=[]
            for x in np.arange(0.001,L+0.1,delta): 
                for y in np.arange(0.001,L+0.1,delta):
                    P=np.array([[x],[y]])
                    f1=self.CRB_RSS_TOA_TDOA_fim(P, RN_RSS, RN_TOA, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TOAStd, TDOAStd)
                    if np.isnan(f1)==0 :
                        CRB.append(f1)

            moy=np.mean(CRB)
            return moy

    def CRB_RSS(self, P, RN_RSS, RSSnp, RSSStd):
            """
            This computes CRB of RSS positioning
            """

            shP     = np.shape(P)
            shRN    = np.shape(RN_RSS)
            RNnum   = shRN[1]

            S = (np.log(10)/10.)* RSSStd/RSSnp

            RNmP = RN_RSS - np.outer(P,np.ones(RNnum))
            mRNmP = (np.sqrt(np.diag(np.dot(RNmP.T,RNmP))))

            num     = np.sum((1+S[:,0]**2)/((S[:,0]**2)*mRNmP**2),axis=0)      # the numerator of the CRLB

            div1    = np.sum(((1+S[:,0]**2)*RNmP**2)/((S[:,0]**2)*mRNmP**4),axis=1).reshape(shP)
            don1    = div1.prod(axis=0)#[0]         # first term of the doniminator

            div2    = np.sum(((1+S[:,0]**2)*(RNmP.prod(axis=0)))/((S[:,0]**2)*mRNmP**4),axis=0)
            don2    = div2**2                       # second term of the doniminator

            CRB     = num/(don1-don2)               # the CRB

            return CRB




    def CRB_TOA(self, P, RN_TOA, TOAStd):
        """
          Compute the CRB in P for the given scenario
        """

        c       = 3e08
        shP     = np.shape(P)
        shRN    = np.shape(RN_TOA)
        RNnum   = shRN[1]
        RoAStd  = c*TOAStd

        num     = np.sum(1/(RoAStd**2),axis=0)[0]  # the numerator of the CRLB

        RNmP    = RN_TOA - outer(P,ones(RNnum))
        mRNmP   = (np.sqrt(diag(dot(RNmP.T,RNmP))))

        div1    = np.sum(RNmP**2/((RoAStd[:,0]**2)*mRNmP**2),axis=1).reshape(shP)
        don1    = div1.prod(axis=0)[0]          # first term of the doniminator

        div2    = np.sum((RNmP.prod(axis=0))/((RoAStd[:,0]**2)*mRNmP**2),axis=0)
        don2    = div2**2                       # second term of the doniminator

        CRB     = num/(don1-don2)               # the CRB

        return CRB

    def CRB_TDOA(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """
            Compute the CRB in P for the given scenario
        """

        c       = 3e08
        shP     = shape(P)
        shRN    = shape(RN1_TDOA)
        RNnum   = shRN[1]
        RDoAStd = c*TDOAStd

        RN1mP   = np.outer(P,np.ones(RNnum))- RN1_TDOA
        mRN1mP  = (np.sqrt(np.diag(np.dot(RN1mP.T,RN1mP)))).reshape(RNnum,1)
        RN2mP   = np.outer(P,np.ones(RNnum))- RN2_TDOA
        mRN2mP  = (np.sqrt(np.diag(np.dot(RN2mP.T,RN2mP)))).reshape(RNnum,1)

        num = np.sum(2/(RDoAStd[:,0]**2)*(1-sum((RN1mP/mRN1mP[:,0])*(RN2mP/mRN2mP[:,0]),axis=0)),axis=0)       # the numerator of the CRLB


        div1 = np.sum((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])**2,axis=1).reshape(shP)
        don1 = div1.prod(axis=0)[0]          # first term of the doniminator

        div2 = np.sum(prod((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0]),axis=0),axis=0)
        don2 = div2**2                       # second term of the doniminator

        CRB = num/(don1-don2)               # the CRB

        return CRB

    def Angle_RSS(self, P, RN_RSS, RSSnp, RSSStd):
        """
        This computes CRB of RSS positioning
        """

        shP     = shape(P)
        shRN    = shape(RN_RSS)
        RNnum   = shRN[1]
        angle   =0.0

        S               = (log(10)/10)* RSSStd/RSSnp

        RNmP    = RN_RSS - outer(P,ones(RNnum))
        mRNmP   = (sqrt(diag(dot(RNmP.T,RNmP))))

        num     = sum((1+S[:,0]**2)/((S[:,0]**2)*mRNmP**2),axis=0)      # the numerator of the CRLB

        for i in range(shRN[1]):
                for j in range(shRN[1]):
                        #angle=angle+((RNmP[0,i]*RNmP[1,j])/(mRNmP[i]*mRNmP[j])-(RNmP[0,j]*RNmP[1,i])/(mRNmP[j]*mRNmP[i]))**2
                        angle=angle+(((RNmP[0,i]*RNmP[1,j])/(mRNmP[i]*mRNmP[j])-(RNmP[0,j]*RNmP[1,i])/(mRNmP[j]*mRNmP[i]))**2)*(((1+S[i,:]**2)*(1+S[j,:]**2))/(2*(S[i,:]*S[j,:]*mRNmP[j]*mRNmP[i])**2))


        return num/angle

    def Angle_TOA(self, P, RN_TOA, TOAStd):
        """
        This computes CRB of RSS positioning
        """
        c       = 3e08
        shP     = shape(P)
        shRN    = shape(RN_TOA)
        RNnum   = shRN[1]
        angle = 0.0
        RNmP    = RN_TOA - outer(P,ones(RNnum))
        mRNmP   = (sqrt(diag(dot(RNmP.T,RNmP))))
        RoAStd  = c*TOAStd
        num     = sum(1/(RoAStd**2),axis=0)[0]
        for i in range(shRN[1]):
                for j in range(shRN[1]):
                        #angle=angle+((RNmP[0,i]*RNmP[1,j])/(mRNmP[i]*mRNmP[j])-(RNmP[0,j]*RNmP[1,i])/(mRNmP[j]*mRNmP[i]))**2
                        angle=angle+(((RNmP[0,i]*RNmP[1,j])/(mRNmP[i]*mRNmP[j])-(RNmP[0,j]*RNmP[1,i])/(mRNmP[j]*mRNmP[i]))**2)/(2*(RoAStd[i,:]*RoAStd[j,:])**2)

        
        
        return num/angle

    def Angle_TDOA(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """
        Compute the CRB in P for the given scenario
        """
        
        c       = 3e08
        shP     = shape(P)
        shRN    = shape(RN1_TDOA)
        RNnum   = shRN[1]
        RDoAStd = c*TDOAStd
        angle=0.0

        RN1mP   = outer(P,ones(RNnum))- RN1_TDOA
        mRN1mP  = (sqrt(diag(dot(RN1mP.T,RN1mP)))).reshape(RNnum,1)
        RN2mP   = outer(P,ones(RNnum))- RN2_TDOA
        mRN2mP  = (sqrt(diag(dot(RN2mP.T,RN2mP)))).reshape(RNnum,1)

        num     = 0.0
        for i in range(shRN[1]):
            num=num+(1-(RN1mP[0,i]*RN2mP[0,i])/(mRN1mP[i]*mRN2mP[i])-(RN1mP[1,i]*RN2mP[1,i])/(mRN1mP[i]*mRN2mP[i]))/RDoAStd[i,0]**2
        div1    = sum((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])**2,axis=1).reshape(shP)
        don1    = div1.prod(axis=0)[0]          # first term of the doniminator

        div2    = sum(prod((1/RDoAStd[:,0]**2)*(-RN1mP/mRN1mP[:,0]+RN2mP/mRN2mP[:,0]),axis=0),axis=0)
        don2    = div2**2

        for i in range(shRN[1]):
            for j in range(shRN[1]):
                anglei1=(RN1mP[0,i]*RN2mP[1,i])/(mRN1mP[i]*mRN2mP[i])-(RN2mP[0,i]*RN1mP[1,i])/(mRN1mP[i]*mRN2mP[i])
                anglej1=(RN1mP[0,j]*RN2mP[1,j])/(mRN1mP[j]*mRN2mP[j])-(RN2mP[0,j]*RN1mP[1,j])/(mRN1mP[j]*mRN2mP[j])
                angleji=(RN1mP[0,i]*RN1mP[1,j])/(mRN1mP[i]*mRN1mP[j])-(RN1mP[0,j]*RN1mP[1,i])/(mRN1mP[i]*mRN1mP[j])
                        
                angle=angle+((angleji-anglej1+anglei1)**2)/(RDoAStd[i,0]*RDoAStd[j,0])**2
        '''print RN1_TDOA
        print RN2_TDOA
        print 2*num/(don1-don2)
        print (2*num)/(angle*0.5)
        print self.CRB_TDOA_fim(P, RN1_TDOA, RN2_TDOA, TDOAStd)'''
        
        return num/angle
            
            
    def P_CRB_RSS(self, P0, RN_RSS, RSSnp, RSSStd):
        """
        This uses RSS CRB to compute the best position of additional reference node
        """

        #P0=array([[0.0],[0.0]])
        P=optimize.fmin(self.CRB_RSS,P0,args=(RN_RSS, RSSnp, RSSStd),xtol=1e-10,ftol=1e-10)


        return P

    def P_CRB_TOA(self, P0, RN_TOA, TOAStd):
        """
        This uses TOA CRB to compute the best position of additional reference node
        """

        #P0=array([[0.0],[0.0]])
        P=optimize.fmin(self.CRB_TOA,P0,args=(RN_TOA, TOAStd),xtol=1e-10,ftol=1e-10)


        return P                


