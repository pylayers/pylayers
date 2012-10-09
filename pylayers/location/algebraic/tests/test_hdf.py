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
#Bernard UGUEN          : buguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
from numpy import *
from scipy import *
from pylayers.location.algebraic.hdf import *
from pylayers.location.algebraic.rss import *
from pylayers.util.PyGraphTool import *
from pylayers.util.PyMathTool import *
import matplotlib.pyplot as plt
import pdb


def plot_cdf_All():
    for i in range(Ntrial):
        print "                           ", Ntrial - i
        P = L * rand(2, 1)
        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
   #     TDoAStd=(sig1/c)*rand(RNnumTDOA,1)
        TDoAStd = (sig1 / c) * ones((RNnumTDOA, 1))
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))

        print RN_RSS, '\n', P, '\n', PL0, '\n', d0, '\n', RSSnp, '\n', RSSStd
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
#                ToAStd=(sig/c)*rand(RNnumTOA,1)
        ToAStd = (sig / c) * ones((RNnumTOA, 1))

        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        print '############### RSS ##################'
        print 'RNRSS\n', RN_RSS
        print    'PL0\n', PL0
        print   'd0\n', d0
        print   'RSS\n', RSS
        print   'RSSnp\n', RSSnp
        print   'RSSSTD\n', RSSStd

        print '############### TOA ##################'
        print 'RNTOA\n', RN_TOA
        print 'ToA\n', ToA
        print   'ToAStd\n', ToAStd

#               print '############### TDOA ##################'
#               print 'RNTDOA\n', RN_TDOA
#               print 'RNTDOA_ref\n', RN_TDOA2
#               print   'TDOA\n', TDoA*1e-9
#               print   'TDOASTD\n', TDoAStd*1e-9

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)

        '''P1=S1.LSHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        LS_vect.append(dist(P1,P))

        P1tls=S1.TLSHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2,
            ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TLS_vect.append(dist(P1tls,P))

        P2=S1.WLSHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        WLS_vect.append(dist(P2,P))'''

        P2twls = S1.TWLSHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
                                  ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLS_vect.append(dist(P2twls, P))

        P0 = L * rand(2, 1)
        P4 = S1.MLHDFLocate(P, P0, RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
                            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        ML_vect.append(dist(P4, P))

        P5 = S1.SDPHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
                             ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDP_vect.append(dist(P5, P))

        '''crb=S1.CRBHDFLocate(P, RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToAStd, TDoAStd, PL0, d0, RSSnp, RSSStd, Rest)
        crb_vect.append(sqrt(crb))'''

    #cdf(LS_vect,"k-","LS",2)
    #cdf(TLS_vect,"k-","TLS",2)
    #cdf(WLS_vect,"r-","WLS",2)

    cdf(ML_vect, "g-", "ML", 2)
    cdf(SDP_vect, "b--", "SDP", 2)
    cdf(TWLS_vect, "r-.", "TWLS", 2)
    #cdf(crb_vect,"g-.","CRLB",2)
    plt.legend(loc=4, numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_RSS_TOA_TDOA.pdf", fromat="pdf")
    plt.close()


def plot_cdf_RSS_TOA():
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)

        '''P1=S1.LSHDFLocate(RN_RSS, RN_TOA, None, None, ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        LS_vect.append(dist(P1,P))'''

        P1tls = S1.TLSHDFLocate(RN_RSS, RN_TOA, None, None, ToA,
                                TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TLS_vect.append(dist(P1tls, P))

        '''P2=S1.WLSHDFLocate(RN_RSS, RN_TOA, None, None, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        WLS_vect.append(dist(P2,P))'''

        '''P2twls=S1.TWLSHDFLocate(RN_RSS, RN_TOA, None, None, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLS_vect.append(dist(P2twls,P))'''

        P0 = L * rand(2, 1)
        P4 = S1.MLHDFLocate(P, P0, RN_RSS, RN_TOA, None, None, ToA,
                            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        ML_vect.append(dist(P4, P))

        P5 = S1.SDPHDFLocate(RN_RSS, RN_TOA, None, None, ToA,
                             ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDP_vect.append(dist(P5, P))

        '''crb=S1.CRBHDFLocate(P, RN_RSS, RN_TOA, None, None, ToAStd, TDoAStd, PL0, d0, RSSnp, RSSStd, Rest)
        crb_vect.append(sqrt(crb))'''

    #cdf(LS_vect,"k-","LS",2)
    #cdf(TLS_vect,"k-","TLS",2)
    #cdf(WLS_vect,"r-","WLS",2)
    cdf(ML_vect, "g-", "ML", 2)
    cdf(SDP_vect, "b--", "SDP", 2)
    cdf(TLS_vect, "r-.", "TWLS", 2)
    #cdf(crb_vect,"g-.","CRLB",2)
    plt.legend(loc=4, numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_RSS_TOA.pdf", fromat="pdf")
    plt.close()


def plot_cdf_RSS_TDOA():
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)

        '''P1=S1.LSHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        LS_vect.append(dist(P1,P))

        P1tls=S1.TLSHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
            TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TLS_vect.append(dist(P1tls,P))

        P2=S1.WLSHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        WLS_vect.append(dist(P2,P))'''

        P2twls = S1.TWLSHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
                                  ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLS_vect.append(dist(P2twls, P))

        P0 = L * rand(2, 1)
        P4 = S1.MLHDFLocate(P, P0, RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
                            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        ML_vect.append(dist(P4, P))

        P5 = S1.SDPHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
                             ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDP_vect.append(dist(P5, P))

        '''crb=S1.CRBHDFLocate(P, RN_RSS, None, RN_TDOA, RN_TDOA2, ToAStd, TDoAStd, PL0, d0, RSSnp, RSSStd, Rest)
        crb_vect.append(sqrt(crb))'''

    #cdf(LS_vect,"k-","LS",2)
    #cdf(TLS_vect,"k-","TLS",2)
    #cdf(WLS_vect,"r-","WLS",2)

    cdf(ML_vect, "g-", "ML", 2)
    cdf(SDP_vect, "b--", "SDP", 2)
    cdf(TWLS_vect, "r-.", "TWLS", 2)
    #cdf(crb_vect,"g-.","CRLB",2)
    plt.legend(loc=4, numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_RSS_TDOA.pdf", fromat="pdf")
    plt.close()


def plot_cdf_TOA_TDOA():
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * ones((RNnumTOA, 1))
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)
        while (sum(ToA < 0)) != 0:
            ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)

        '''P1=S1.LSHDFLocate(None, RN_TOA, RN_TDOA, RN_TDOA2, ToA, TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        LS_vect.append(dist(P1,P))

        P1tls=S1.TLSHDFLocate(None, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
            TDoA, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TLS_vect.append(dist(P1tls,P))

        P2=S1.WLSHDFLocate(None, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        WLS_vect.append(dist(P2,P))'''
        print '############### TOA ##################'
        print 'RNTOA\n', RN_TOA
        print 'ToA\n', ToA
        print   'ToAStd\n', ToAStd

#               print '############### TDOA ##################'
#               print 'RNTDOA\n', RN_TDOA
#               print 'RNTDOA_ref\n', RN_TDOA2
#               print   'TDOA\n', TDoA
#               print   'TDOASTD\n', TDoAStd
#                P2twls=S1.TWLSHDFLocate(None, RN_TOA, RN_TDOA, RN_TDOA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
#                TWLS_vect.append(dist(P2twls,P))
#               ToAStdsave.append(ToAStd)
#               print 'point:',P2twls

        P0 = L * rand(2, 1)
        P4 = S1.MLHDFLocate(P, P0, None, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
                            ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)

        ML_vect.append(dist(P4, P))

#                P5=S1.SDPHDFLocate(None, RN_TOA, RN_TDOA, RN_TDOA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
#                SDP_vect.append(dist(P5,P))

        '''crb=S1.CRBHDFLocate(P, None, RN_TOA, RN_TDOA, RN_TDOA2, ToAStd, TDoAStd, PL0, d0, RSSnp, RSSStd, Rest)
        crb_vect.append(sqrt(crb))'''

    #cdf(LS_vect,"k-","LS",2)
    #cdf(TLS_vect,"k-","TLS",2)
    #cdf(WLS_vect,"r-","WLS",2)

    cdf(ML_vect, "g-", "ML", 2)
    #cdf(SDP_vect,"b--","SDP",2)
    #cdf(TWLS_vect,"r-.","TWLS",2)
    #cdf(crb_vect,"g-.","CRLB",2)
    #plt.legend(loc=4,numpoints=1)
    #plt.axis([0,10,0,1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_TOA_TDOA.pdf", fromat="pdf")
    plt.close()


def plot_cdf_TOA():
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * ones((RNnumTOA, 1))
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        while (sum(ToA < 0)) != 0:
            ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        print '############### TOA ##################'
        print 'RNTOA\n', RN_TOA
        print 'ToA\n', ToA
        print   'ToAStd\n', ToAStd

        P0 = L * rand(2, 1)

        S3 = ToALocation(RN_TOA)
        P4 = S3.MLToALocate(P, P0, RN_TOA, ToA, ToAStd)
        print dist(P4, P)
        ML_vect.append(dist(P4, P))

    cdf(ML_vect, "g-", "ML", 2)
#        #cdf(SDP_vect,"b--","SDP",2)
#        #cdf(TWLS_vect,"r-.","TWLS",2)
#        #cdf(crb_vect,"g-.","CRLB",2)
#        #plt.legend(loc=4,numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.show()
#        plt.savefig("HDF_cdf_TOA_TDOA.pdf", fromat="pdf")
#        plt.close()


def plot_cdf_ML_All():
    MLrss_vect = []
    MLrsstoa_vect = []
    MLrsstdoa_vect = []
    MLrsstoatdoa_vect = []
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        P0 = L * rand(2, 1)
        P4rss = S2.MLDRSSLocate(P, P0, RN_RSS, PL0, d0, RSS, RSSnp, RSSStd)
        MLrss_vect.append(dist(P4rss, P))

        P4rsstoa = S1.MLHDFLocate(P, P0, RN_RSS, RN_TOA, None, None, ToA,
                                  ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrsstoa_vect.append(dist(P4rsstoa, P))

        P4rsstdoa = S1.MLHDFLocate(P, P0, RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
                                   ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrsstdoa_vect.append(dist(P4rsstdoa, P))

        P4rsstoatdoa = S1.MLHDFLocate(P, P0, RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA, ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrsstoatdoa_vect.append(dist(P4rsstoatdoa, P))

    cdf(MLrss_vect, "k-", "RSSI", 2)
    cdf(MLrsstoa_vect, "b-", "RSSI + TOA", 2)
    cdf(MLrsstdoa_vect, "r-", "RSSI + TDOA", 2)
    cdf(MLrsstoatdoa_vect, "g-", "RSSI + TOA + TDOA", 2)

    plt.legend(loc=4, numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_ML.pdf", fromat="pdf")
    plt.close()


def plot_cdf_TWLS_All():
    TWLSrss_vect = []
    TWLSrsstoa_vect = []
    TWLSrsstdoa_vect = []
    TWLSrsstoatdoa_vect = []
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        P4rss = S2.TWLSRSSLocate(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLSrss_vect.append(dist(P4rss, P))

        P4rsstoa = S1.TWLSHDFLocate(RN_RSS, RN_TOA, None, None, ToA,
                                    ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLSrsstoa_vect.append(dist(P4rsstoa, P))

        P4rsstdoa = S1.TWLSHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
                                     ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLSrsstdoa_vect.append(dist(P4rsstdoa, P))

        P4rsstoatdoa = S1.TWLSHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
                                        ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        TWLSrsstoatdoa_vect.append(dist(P4rsstoatdoa, P))

    cdf(TWLSrss_vect, "k-", "RSSI", 2)
    cdf(TWLSrsstoa_vect, "b-", "RSSI + TOA", 2)
    cdf(TWLSrsstdoa_vect, "r-", "RSSI + TDOA", 2)
    cdf(TWLSrsstoatdoa_vect, "g-", "RSSI + TOA + TDOA", 2)

    plt.legend(loc=4, numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_TWLS.pdf", fromat="pdf")
    plt.close()


def plot_cdf_SDP_All():
    SDPrss_vect = []
    SDPrsstoa_vect = []
    SDPrsstdoa_vect = []
    SDPrsstoatdoa_vect = []
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        P4rss = S2.SDPRSSLocate(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss_vect.append(dist(P4rss, P))

        P4rsstoa = S1.SDPHDFLocate(RN_RSS, RN_TOA, None, None, ToA,
                                   ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrsstoa_vect.append(dist(P4rsstoa, P))

        P4rsstdoa = S1.SDPHDFLocate(RN_RSS, None, RN_TDOA, RN_TDOA2, ToA,
                                    ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrsstdoa_vect.append(dist(P4rsstdoa, P))

        P4rsstoatdoa = S1.SDPHDFLocate(RN_RSS, RN_TOA, RN_TDOA, RN_TDOA2, ToA,
                                       ToAStd, TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrsstoatdoa_vect.append(dist(P4rsstoatdoa, P))

    cdf(SDPrss_vect, "k-", "RSSI", 2)
    cdf(SDPrsstoa_vect, "b-", "RSSI + TOA", 2)
    cdf(SDPrsstdoa_vect, "r-", "RSSI + TDOA", 2)
    cdf(SDPrsstoatdoa_vect, "g-", "RSSI + TOA + TDOA", 2)

    plt.legend(loc=4, numpoints=1)
    plt.axis([0, 10, 0, 1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_SDP.pdf", fromat="pdf")
    plt.close()


def plot_cdf_SDP_RSS_nTOA():
    SDPrss_vect = []
    SDPrss1toa_vect = []
    SDPrss2toa_vect = []
    SDPrss3toa_vect = []
    SDPrss4toa_vect = []
    #Ntrial=10
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        P4rss = S2.SDPRSSLocate(RN_RSS, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss_vect.append(dist(P4rss, P))

        P4rss1toa = S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:, 0:1]), None, None, ToA[0:1, :], ToAStd[0:1, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss1toa_vect.append(dist(P4rss1toa, P))

        P4rss2toa = S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:, 0:2]), None, None, ToA[0:2, :], ToAStd[0:2, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss2toa_vect.append(dist(P4rss2toa, P))

        P4rss3toa = S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:, 0:3]), None, None, ToA[0:3, :], ToAStd[0:3, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss3toa_vect.append(dist(P4rss3toa, P))

        P4rss4toa = S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:, 0:4]), None, None, ToA[0:4, :], ToAStd[0:4, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss4toa_vect.append(dist(P4rss4toa, P))

    cdf(SDPrss_vect, "k-", "RSSI", 2)
    cdf(SDPrss1toa_vect, "b-", "RSSI + 1 TOA", 2)
    cdf(SDPrss2toa_vect, "r-", "RSSI + 2 TOA", 2)
    cdf(SDPrss3toa_vect, "y-", "RSSI + 3 TOA", 2)
    cdf(SDPrss4toa_vect, "g-", "RSSI + 4 TOA", 2)

    plt.legend(loc=4, numpoints=1)
    #plt.axis([0,10,0,1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_SDP_RSS_nTOA.pdf", fromat="pdf")
    plt.close()

    ntoa = range(5)
    err_moy = [mean(SDPrss_vect), mean(SDPrss1toa_vect), mean(
        SDPrss2toa_vect), mean(SDPrss3toa_vect), mean(SDPrss4toa_vect)]
    plt.plot(ntoa, err_moy, "ko-", linewidth=2)
    plt.xlim = 5
    plt.grid('on')
    plt.xlabel("Number of Added TOA")
    plt.ylabel("Average Positioning Error")
    plt.savefig("HDF_err_SDP_RSS_nTOA.pdf", fromat="pdf")
    plt.close()


def plot_cdf_SDP_RSS_nTOA_crlb():
    SDPrss_vect = []
    SDPrss1toa_vect = []
    SDPrss1toacrlb_vect = []
    SDPrss2toa_vect = []
    SDPrss3toa_vect = []
    SDPrss4toa_vect = []
    #Ntrial=10
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        S3 = CRBLocation(RN_RSS)
        CRB_RSS_TOA_fim = []
        a = (4 * rand(1)).astype(int)[0]
        for b in range(3):
            CRB_RSS_TOA_fim.append(sqrt(S3.CRB_RSS_TOA_fim(P, RN_RSS, RN_TOA[
                :, b:b + 1], RSSnp, RSSStd, ToAStd[b:b + 1, :])))
        u = where(CRB_RSS_TOA_fim == min(CRB_RSS_TOA_fim))
        v = where(CRB_RSS_TOA_fim == max(CRB_RSS_TOA_fim))
        b = u[0][0]
        a = v[0][0]
        P4rss1toa = S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:, a:a + 1]), None, None, ToA[a:a + 1, :], ToAStd[a:a + 1, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss1toa_vect.append(dist(P4rss1toa, P))

        P4rss1toacrlb = S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:, b:b + 1]), None, None, ToA[b:b + 1, :], ToAStd[b:b + 1, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss1toacrlb_vect.append(dist(P4rss1toacrlb, P))

        '''P4rss2toa=S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:,0:2]), None, None, ToA[:,0:2], ToAStd[:,0:2], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss2toa_vect.append(dist(P4rss2toa,P))

        P4rss3toa=S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:,0:3]), None, None, ToA[:,0:3], ToAStd[:,0:3], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss3toa_vect.append(dist(P4rss3toa,P))

        P4rss4toa=S1.SDPHDFLocate(RN_RSS, matrix(RN_TOA[:,0:4]), None, None, ToA[:,0:4], ToAStd[:,0:4], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        SDPrss4toa_vect.append(dist(P4rss4toa,P))'''

    cdf(SDPrss1toa_vect, "b-", "RSSI + 1 TOA (Randomly)", 2)
    cdf(SDPrss1toacrlb_vect, "r-", "RSSI + 1 TOA (CRLB Criteria)", 2)

    plt.legend(loc=4, numpoints=1)
    #plt.axis([0,10,0,1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_SDP_RSS_nTOA_crlb.pdf", fromat="pdf")
    plt.close()


def plot_cdf_ML_RSS_nTOA():
    MLrss_vect = []
    MLrss1toa_vect = []
    MLrss2toa_vect = []
    MLrss3toa_vect = []
    MLrss4toa_vect = []
    #Ntrial=10
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        P0 = L * rand(2, 1)
        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        P4rss = S2.MLDRSSLocate(P, P0, RN_RSS, PL0, d0, RSS, RSSnp, RSSStd)
        MLrss_vect.append(dist(P4rss, P))

        P4rss1toa = S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:, 0:1]), None, None, ToA[0:1, :], ToAStd[0:1, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss1toa_vect.append(dist(P4rss1toa, P))

        P4rss2toa = S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:, 0:2]), None, None, ToA[0:2, :], ToAStd[0:2, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss2toa_vect.append(dist(P4rss2toa, P))

        P4rss3toa = S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:, 0:3]), None, None, ToA[0:3, :], ToAStd[0:3, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss3toa_vect.append(dist(P4rss3toa, P))

        P4rss4toa = S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:, 0:4]), None, None, ToA[0:4, :], ToAStd[0:4, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss4toa_vect.append(dist(P4rss4toa, P))

    cdf(MLrss_vect, "k-", "RSSI", 2)
    cdf(MLrss1toa_vect, "b-", "RSSI + 1 TOA", 2)
    cdf(MLrss2toa_vect, "r-", "RSSI + 2 TOA", 2)
    cdf(MLrss3toa_vect, "y-", "RSSI + 3 TOA", 2)
    cdf(MLrss4toa_vect, "g-", "RSSI + 4 TOA", 2)

    plt.legend(loc=4, numpoints=1)
    #plt.axis([0,10,0,1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_ML_RSS_nTOA.pdf", fromat="pdf")
    plt.close()

    ntoa = range(5)
    err_moy = [mean(MLrss_vect), mean(MLrss1toa_vect), mean(
        MLrss2toa_vect), mean(MLrss3toa_vect), mean(MLrss4toa_vect)]
    plt.plot(ntoa, err_moy, "ko-", linewidth=2)
    plt.xlim = 5
    plt.grid('on')
    plt.xlabel("Number of Added TOA")
    plt.ylabel("Average Positioning Error")
    plt.savefig("HDF_err_ML_RSS_nTOA.pdf", fromat="pdf")
    plt.close()


def plot_cdf_ML_RSS_nTOA_crlb():
    MLrss_vect = []
    MLrss1toa_vect = []
    MLrss1toacrlb_vect = []
    MLrss2toa_vect = []
    MLrss3toa_vect = []
    MLrss4toa_vect = []
    #Ntrial=10
    for i in range(Ntrial):
        print "                                                                            ", Ntrial - i
        P = L * rand(2, 1)

        shRN_TDOA = shape(RN_TDOA)
        RNnumTDOA = shRN_TDOA[1]
        RNTDOAmp = RN_TDOA - P
        RNTDOA2mp = RN_TDOA2 - P
        RNTDOAmp2 = (
            sum(RNTDOAmp * RNTDOAmp, axis=0)).reshape(RNnumTDOA, 1)
        RNTDOA2mp2 = (
            sum(RNTDOA2mp * RNTDOA2mp, axis=0)).reshape(RNnumTDOA, 1)
        RDoA = sqrt(RNTDOAmp2) - sqrt(RNTDOA2mp2)
        TDoAStd = (sig1 / c) * rand(RNnumTDOA, 1)
        TDoA = RDoA / c + TDoAStd * randn(RNnumTDOA, 1)

        shRN_RSS = shape(RN_RSS)
        RNnumRSS = shRN_RSS[1]
        RSSStd = sh * ones((RNnumRSS, 1))
        RSSnp = np * ones((RNnumRSS, 1))
        d0 = 1.0
        S2 = RSSLocation(RN_RSS)
        PL0 = pl0 * ones((RNnumRSS, 1))
        RSS = S2.getPL(RN_RSS, P, PL0, d0, RSSnp, RSSStd)

        shRN_TOA = shape(RN_TOA)
        RNnumTOA = shRN_TOA[1]
        RNTOAmp = RN_TOA - P
        RNTOAmp2 = (sum(RNTOAmp * RNTOAmp, axis=0)).reshape(RNnumTOA, 1)
        RoA = sqrt(RNTOAmp2)
        ToAStd = (sig / c) * rand(RNnumTOA, 1)
        ToA = RoA / c + ToAStd * randn(RNnumTOA, 1)

        P0 = L * rand(2, 1)
        S1 = HDFLocation(RN_RSS, RN_TOA, RN_TDOA)
        S2 = RSSLocation(RN_RSS)
        S3 = CRBLocation(RN_RSS)
        CRB_RSS_TOA_fim = []
        a = (4 * rand(1)).astype(int)[0]
        for b in range(4):
            CRB_RSS_TOA_fim.append(sqrt(S3.CRB_RSS_TOA_fim(P, RN_RSS, RN_TOA[
                :, b:b + 1], RSSnp, RSSStd, ToAStd[b:b + 1, :])))
        u = where(CRB_RSS_TOA_fim == min(CRB_RSS_TOA_fim))
        v = where(CRB_RSS_TOA_fim == max(CRB_RSS_TOA_fim))
        b = u[0][0]
        a = v[0][0]
        P4rss1toa = S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:, a:a + 1]), None, None, ToA[a:a + 1, :], ToAStd[a:a + 1, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss1toa_vect.append(dist(P4rss1toa, P))

        P4rss1toacrlb = S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:, b:b + 1]), None, None, ToA[b:b + 1, :], ToAStd[b:b + 1, :], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss1toacrlb_vect.append(dist(P4rss1toacrlb, P))

        '''P4rss2toa=S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:,0:2]), None, None, ToA[:,0:2], ToAStd[:,0:2], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss2toa_vect.append(dist(P4rss2toa,P))

        P4rss3toa=S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:,0:3]), None, None, ToA[:,0:3], ToAStd[:,0:3], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss3toa_vect.append(dist(P4rss3toa,P))

        P4rss4toa=S1.MLHDFLocate(P, P0, RN_RSS, matrix(RN_TOA[:,0:4]), None, None, ToA[:,0:4], ToAStd[:,0:4], TDoA, TDoAStd, PL0, d0, RSS, RSSnp, RSSStd, Rest)
        MLrss4toa_vect.append(dist(P4rss4toa,P))'''

    cdf(MLrss1toa_vect, "b-", "RSSI + 1 TOA (Randomly)", 2)
    cdf(MLrss1toacrlb_vect, "r-", "RSSI + 1 TOA (CRLB Criteria)", 2)

    plt.legend(loc=4, numpoints=1)
    #plt.axis([0,10,0,1])
    plt.grid('on')
    plt.xlabel("Positioning error (m)")
    plt.ylabel("Cumulative probability")
    plt.savefig("HDF_cdf_ML_RSS_nTOA_crlb.pdf", fromat="pdf")
    plt.close()


if (__name__ == "__main__"):

    Ntrial = 500
    LS_vect = []
    TLS_vect = []
    WLS_vect = []
    TWLS_vect = []
    TS_vect = []
    ML_vect = []
    SDP_vect = []
    crb_vect = []
    sig_vect = []

    Rest = 'mode'
    c = 3e08
    L = 20.0

    '''RN_RSS=array([[0.0,L,0.0,L],[0.0,0.0,L,L]])
    RN_TOA=array([[L/2.0,3*L/4.0,L/2.0,L/4.0],[L/4.0,L/2,3*L/4.0,L/2.0]])
    RN_TDOA_all=array([[L/2.0,2*L,L/2.0,-L],[-L,L/2.0,2*L,L/2.0]])
    RN_TDOA=array([[2*L,L/2.0,-L],[L/2.0,2*L,L/2.0]])
    RN_TDOA2=RN_TDOA_all[:,0:1]*ones(shape(RN_TDOA))'''

    RN_RSS = array([[0.0, L, 0.0, L], [0.0, 0.0, L, L]])
    RN_TOA = array(
        [[L / 3.0, L, 2 * L / 3.0, 0.0], [0.0, L / 3.0, L, 2 * L / 3.0]])
    RN_TDOA_all = array(
        [[L / 2.0, 2 * L, L / 2.0, -L], [-L, L / 2.0, 2 * L, L / 2.0]])
    RN_TDOA = array([[2 * L, L / 2.0, -L], [L / 2.0, 2 * L, L / 2.0]])
    RN_TDOA2 = RN_TDOA_all[:, 0:1] * ones(shape(RN_TDOA))

    RN = [RN_RSS.T, RN_TOA.T, RN_TDOA_all.T]

    '''plt.plot(RN_RSS[0,:],RN_RSS[1,:],'ro')
    plt.plot(RN_TOA[0,:],RN_TOA[1,:],'go')
    plt.plot(RN_TDOA_all[0,:],RN_TDOA_all[1,:],'bo')
    plt.show()'''

    pl0 = -34.7
    np = 2.645
    sh = 4.34
    sig = 2.97
    sig1 = 3.55
    ToAStdsave = []
    plot_cdf_All()
#    plot_cdf_RSS_TOA()
#    plot_cdf_RSS_TDOA()
#    plot_cdf_TOA_TDOA()
#    plot_cdf_SDP_All()
#    plot_cdf_ML_All()
#    plot_cdf_ML_RSS_nTOA()
#    plot_cdf_ML_RSS_nTOA_crlb()
#    plot_cdf_SDP_RSS_nTOA()
#    plot_cdf_SDP_RSS_nTOA_crlb()
