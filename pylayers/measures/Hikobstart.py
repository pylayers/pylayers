# -*- coding:Utf-8 -*-
import os
import pdb
import sys
import pandas as pd
import numpy as np
import scipy.io as io
from pylayers.util.project import *
import matplotlib.pyplot as plt
import copy

class Hikob(PyLayers):
    """
    """
    def load(day=11,serie='',scenario='20',run=1):
        if day==11:
            files = os.listir(CORMORAN+'/POST-TREATED/11-06-2014/HIKOB')
        if day==12:
            files = os.listdir(CORMORAN+'/POST-TREATED/12-06-2014/HIKOB')

        if serie != '':
            filename = filter(lambda x : 'S'+serie in x ,files)
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            filename = filter(lambda x : 'R'+str(run) in x ,filsc)

        data = io.loadmat(filename)
        self.rssi = data['rssi']
        self.t = data['ts']



def extract(D):
    # conversion pandas -> numpy
    seq = np.array(D['seq'])
    rssi = np.array(D['rssi'])
    #
    # find start of sequence seq==1
    #
    u = np.where(seq==1)[0]
    trssi = {}
    for k in range(len(u)):
        if k!=(len(u)-1):
            trssi[k]=rssi[u[k]:u[k+1]]
        else:
            trssi[k]=rssi[u[k]:]
    return(trssi)

CORMORAN='/home/uguen/svn2/measures/CORMORAN'
os.chdir(CORMORAN+'/RAW/11-06-2014/HIKOB')
files = os.listdir('.')
rawfiles = filter(lambda x: 'RAW' in x ,files)

HKB ={'AP1':1,'AP2':2,'AP3':3,'AP4':4,'HeadRight':5,'TorsoTopRight':6,'TorsoTopLeft':7,'BackCenter':8,'ElbowRight':9,'ElbowLeft':10,'HipRight':11,'WristRight':12,'WristLeft':13,'KneeLeft':14,'AnckleRight':16,'AnckleLeft':15}
for k in HKB:
    print k,HKB[k]
JiHadHKB = {'TorsoTopRight':10,'TorsoTopLeft':9,'BackCenter':11,'ShoulderLeft':12}
NicolasHKB = {'TorsoTopRight':6,'TorsoTopLeft':5,'BackCenter':7,'ShoulderLeft':8}
EricKKB = {'TorsoTopRight':15,'TorsoTopLeft':13,'BackCenter':16,'ShoulderLeft':14}

links = [176,75,76,77,78]


data = {}

dt1=[('src','uint8'),
    ('seq','uint32'),
    ('rssi','int8'),
    ('lqi','uint8'),
    ('n_sensors','uint8'),
    ('data_length','uint8'),
    ('interslot','uint16'),
    ('other','int8')]

rawfiles.sort()
#pdb.set_trace()
for k,fi in enumerate(rawfiles):
    dt = np.fromfile(fi,dtype=dt1)
    D = pd.DataFrame()
    D = D.from_records(dt)
    del D['lqi']
    del D['data_length']
    del D['other']
    del D['n_sensors']
    del D['interslot']
    data[k+1]=D

# data extraction
for k in range(1,17):
    for l in range(1,17):
        dt = data[k][data[k]['src']==l]['seq']
        udt = np.where(dt<2)

        Nsc = len(np.where(np.diff(udt)!=1)[0])

data2 = {}
for k in range(1,17):
    data2[k] = {}
    for l in range(1,17):
        D=data[k][data[k]['src']==l]
        trssi = extract(D)
        data2[k][l] = trssi

#
# For diagonal terms the rssi table is set to -60
#
tA = {}
A  = []


#
# Estimate median time duration of each series
#
scduration = {}
for sc in range(25):
    scduration[sc] = []
    for k in range(1,17):
        for l in range(1,17):
            if k!=l:
                try:
                    scduration[sc].append(len(data2[k][l][sc]))
                except:
                    pass
            else:
                scduration[sc].append(0)
dti = {}
for k in scduration:
    dti[k]=np.median(np.array(scduration[k])).astype('int')



tstart=np.zeros((16,16))
data2 ={}
for k in range(1,17):
    data2[k]={}
    for l in range(1,17):
        if k!=l:
            dt = data[k][data[k]['src']==l]
            seq  = dt['seq'].values
            rssi =  dt['rssi'].values
            u1 = np.where(seq==1)[0][1]
            u2 = np.where(seq==1)[0][2]
            if u2 > 1000:
                u2 = u1
            tstart[k-1][l-1] = u2
            data2[k][l] = dt[u2:]
        else:
            pass

def matstart(data,index):
    tstart=np.zeros((16,16))
    for k in range(1,17):
        for l in range(1,17):
            if k!=l:
                dt = data[k][l]
                seq  = dt['seq'].values
                rssi =  dt['rssi'].values
                u = np.where(seq==1)[0][index]
                tstart[k-1][l-1] = u.astype('int')
    return tstart


def extracts(data2,tstart):
    datao ={}
    serie ={}
    data2c = copy.deepcopy(data2)
    for k in range(1,17):
        datao[k]={}
        for l in range(1,17):
            if k!=l:
                dt = data2c[k][l]
                datao[k][l] = copy.deepcopy(dt[tstart[k-1][l-1]:])
                dt.drop_duplicates('seq',inplace=True)
                rssi = dt[0:tstart[k-1][l-1]]['rssi'].values
                ind  = dt[0:tstart[k-1][l-1]]['seq'].values

                serie[(k,l)] = pd.Series(rssi,index=ind)
                #serie[k][l].index = serie[k][l]['seq']
                #del serie[k][l]['src']
                #del serie[k][l]['seq']
    return(serie,datao)


def plotindex(dd):
    f,ax=plt.subplots(nrows=16,ncols=16)
    
    for k in range(1,17):
        for l in range(1,17):
            if k!=l:
                ax[k-1,l-1].plot(dd[k,l].index)

def plotindex2(dd):
    f,ax=plt.subplots(nrows=1,ncols=1)
    
    for k in range(1,17):
        for l in range(1,17):
            if k!=l:
                ax.plot(dd[k,l].index,'.')

def plotdata(dd):
    f,ax=plt.subplots(nrows=16,ncols=16)
    
    for k in range(1,17):
        for l in range(1,17):
            if k!=l:
                dd[k][l].plot(ax=ax[k-1,l-1])

def check(dd,ss,ii=1):
    X={}
    for k in range(1,17):
        X[k]={}
        for l in range(1,17):
            if k!=l:
                X[k][l]=np.where(dd[k][l][0:ss[k-1,l-1]]['seq']==ii)[0]
    return X

def dupli(dd,ss):
    X={}
    for k in range(1,17):
        X[k]={}
        for l in range(1,17):
            if k!=l:
                X[k][l]=np.where(dd[k][l][0:ss[k-1,l-1]]['seq'].duplicated())
    return X

def drop_dupli(df):
    for k in range(1,17):
        for l in range(1,17):
            if k!=l:
                df[k,l].drop_duplicates('seq',inplace=True)

def toarray(S):
    """ convert data frame to array
    """
    U  = np.array(S)
    Nt = U.shape[0]
    A  = np.empty((Nt,256))
    A[:] = np.NaN
    u  = np.arange(256)
    v  = np.arange(0,256,17)
    w  = np.setdiff1d(u,v)
    A[:,w] = U
    A = A.reshape(Nt,16,16)
    return(A)

print "get serie 2"

tstart = matstart(data2,1).astype('int')
tstart[6,8] = 3678
tstart[8,6] = 3701

s,data3 = extracts(data2,tstart)
S2 = pd.DataFrame(s)


print "get serie 3"

tstart3 = matstart(data3,1).astype('int')
tstart3[6,8] = 3647
tstart3[8,6] = 3663

s,data4 = extracts(data3,tstart3)

S3 = pd.DataFrame(s)

print "get serie 4 --"

tstart4 = matstart(data4,1).astype('int')
tstart4[9,10] = 93
tstart4[14,8] = 75

s,data5 = extracts(data4,tstart4)

S4 = pd.DataFrame(s)


print "get serie 5 --"

tstart5 = matstart(data5,1).astype('int')
tstart5[9,10] = 48
tstart5[14,8] = 22

s,data6 = extracts(data5,tstart5)

S5 = pd.DataFrame(s)

print "get serie 6 --"

tstart6 = matstart(data6,1).astype('int')
tstart6[9,10] = 48
tstart6[14,8] = 45
tstart6[8,14] = 48

s,data7 = extracts(data6,tstart6)

S6 = pd.DataFrame(s)


print "get serie 7 --"

tstart7 = matstart(data7,1).astype('int')
tstart7[9,10] = 37
tstart7[14,8] = 27
tstart7[8,14] = 35

s,data8 = extracts(data7,tstart7)

S7 = pd.DataFrame(s)
##########
print "get serie 8"

tstart8 = matstart(data8,1).astype('int')
s,data9 = extracts(data8,tstart8)

S8 = pd.DataFrame(s)

print "get serie 9"

tstart9 = matstart(data9,1).astype('int')
s,data10 = extracts(data9,tstart9)
S9 = pd.DataFrame(s)

print "get serie 10"

tstart10 = matstart(data10,1).astype('int')
#plt.plot(s[11,7].index)
tstart10[8,6] = 4316
tstart10[10,6] = 4400
s,data11 = extracts(data10,tstart10)
#plotindex(s)
S10 = pd.DataFrame(s)


print "get serie 11"

tstart11 = matstart(data11,1).astype('int')
tstart11[1,6] = 4332
tstart11[8,6] = 4209
tstart11[10,6] = 4289
s,data12 = extracts(data11,tstart11)

S11 = pd.DataFrame(s)


print "get serie 12"

tstart12 = matstart(data12,1).astype('int')
tstart12[1,6] = 4110
s,data13 = extracts(data12,tstart12)
S12 = pd.DataFrame(s)

print "get serie 13"

tstart13 = matstart(data13,1).astype('int')
s,data14 = extracts(data13,tstart13)
S13 = pd.DataFrame(s)


print "get serie 14"

tstart14 = matstart(data14,1).astype('int')
tstart14[7,5]=3030
tstart14[9,5]=2974
tstart14[9,6]=3136
s,data15 = extracts(data14,tstart14)
S14 = pd.DataFrame(s)


print "get serie 15"

tstart15 = matstart(data15,1).astype('int')
tstart15[7,5]=4036
tstart15[9,5]=4013
tstart15[9,6]=4151

s,data16 = extracts(data15,tstart15)
S15 = pd.DataFrame(s)

print "get serie 16"

tstart16 = matstart(data16,1).astype('int')
s,data17 = extracts(data16,tstart16)
S16 = pd.DataFrame(s)

print "get serie 17"

tstart17 = matstart(data17,1).astype('int')
s,data18 = extracts(data17,tstart17)
S17 = pd.DataFrame(s)

print "get serie 18"

tstart18 = matstart(data18,1).astype('int')
s,data19 = extracts(data18,tstart18)
S18 = pd.DataFrame(s)

print "get serie 19"

tstart19 = matstart(data19,1).astype('int')
tstart19[0,5]=5112
tstart19[3,5]=5104
tstart19[4,5]=5113
tstart19[6,8]=5027
tstart19[6,5]=5109
tstart19[7,5]=5008
tstart19[8,5]=5098
tstart19[8,6]=5035
tstart19[9,5]=4929
tstart19[10,5]=5114
tstart19[12,5]=5092
tstart19[13,5]=5069
tstart19[14,5]=5102
tstart19[15,5]=5108
s,data20 = extracts(data19,tstart19)
S19 = pd.DataFrame(s)


def findmin(ss,value,med):
    uss = np.array(np.where(ss>value))
    for k in range(len(uss)):
        mss = np.min(ss[uss[k][0],uss[k][1]][med-100:med+100])
        ss[uss[k][0],uss[k][1]]= np.where(ss[uss[k][0],uss[k][1]]==mss)[0]
    return ss

print "get serie 20"

tstart20 = matstart(data20,1).astype('int')
tstart20[0,5] = 3972
tstart20[3,5] = 3964
tstart20[4,5] = 3978
tstart20[6,5] = 3975
tstart20[6,8] = 3958
tstart20[6,11] =3580
tstart20[7,5] = 3728
tstart20[8,5] = 3968
tstart20[8,6] = 3858
tstart20[9,5] = 3750
tstart20[10,5] = 3978
tstart20[11,6] = 3574
tstart20[12,5] = 3970
tstart20[13,5] = 3928
tstart20[14,4] = 3570
tstart20[14,5] = 3964
tstart20[15,5] = 3964
s,data21 = extracts(data20,tstart20)
S20 = pd.DataFrame(s)

print "get serie 21"
tstart21 = matstart(data21,1).astype('int')
tstart21[6,8]=7808
tstart21[6,11]=4820
tstart21[11,6]=4817
tstart21[14,4]=4751

s,data22=extracts(data21,tstart21)
S21 = pd.DataFrame(s)

print "get serie 22"

tstart22 = matstart(data22,1).astype('int')
s,data22=extracts(data22,tstart22)
S22 = pd.DataFrame(s)

print "get serie 23"

tstart23 = matstart(data22,1).astype('int')
s,data23=extracts(data23,tstart23)
S23 = pd.DataFrame(s)
#
#
#for sc in range(25):
#    del A
#    cpt = 0
#    for k in range(1,17):
#        for l in range(1,17):
#            #if k!=l:
#            cpt = cpt+1
#            try:
#                rssi = data2[k][l][sc]
#            except:
#                rssi = -100*np.ones(dti[sc])
#            if cpt==1:
#                dim = len(rssi)
#            try:
#                if (len(rssi) < dim):
#                    u = np.zeros(dim)
#                    u[0:len(rssi)]=rssi
#                    rssi = u
#                if (len(rssi) > dim):
#                    rssi = rssi[0:dim]
#                A = np.vstack((A,rssi))
#            except:
#                A = rssi
#
#
#    tA[sc]=A
#
#
#iHKB={}
#for k in HKB:
#    iHKB[HKB[k]]=k
#
#
#dHKB={}
#cpt = 0
#for k in range(1,17):
#    for l in range(1,17):
#        #if k!=l:
#        #dHKB[cpt]=iHKB[k]+' - '+iHKB[l]
#            dHKB[(k,l)]=iHKB[k]+' - '+iHKB[l]
#            cpt = cpt + 1
#
#
#for sc in tA.keys():
#    print sc,tA[sc].shape[1]*25.832e-3/60.
#
#
#dHKBc={}
#dHKBc[2]='Sc20_S5_R1_HKBS_8'
#dHKBc[3]='Sc20_S6_R2_HKBS_9'
#dHKBc[8]='Sc21a_S13_R1_HKBS_16'
#dHKBc[9]='Sc21a_S14_R2_HKBS_17'
#dHKBc[10]='Sc21a_S15_R3_HKBS_18'
#dHKBc[11]='Sc21a_S16_R4_HKBS_19'
#dHKBc[12]='Sc21b_S21_R1_HKBS_24'
#dHKBc[13]='Sc21b_S22_R2_HKBS_25'
#dHKBc[14]='Sc21b_S23_R3_HKBS_26'
#dHKBc[15]='Sc21b_S24_R4_HKBS_27'
#dHKBc[16]='Sc21d_S27_R1_HKBS_31'
#dHKBc[17]='Sc21d_S28_R2_HKBS_33'
#dHKBc[18]='Sc21e_S29_R1_HKBS_34'
#dHKBc[19]='Sc21e_S30_R2_HKBS_35'
#dHKBc[20]='Sc20_S31_R5_HKBS_36'
#dHKBc[21]='Sc21a_S32_R1_Full_38'
#dHKBc[22]='Sc21a_S33_R2_Full_39'
#dHKBc[23]='Sc21a_S34_R3_Full_40'
#dHKBc[24]='Sc21a_S35_R4_Full_41'
#
#
#for k in dHKBc:
#    d = {}
#    Nt = np.shape(tA[k])[1]
#    d['rssi'] = tA[k].reshape(16,16,Nt)
#    d['t'] = np.arange(Nt)*25.832e-3
#    filename = dHKBc[k]+'.mat'
#    io.savemat(filename,d)
#
#
### Sc20_S5_R1_HKBS_8
##links = [176,75,76,77,78]
##plt.figure(figsize=(20,5))
##
##for l in links:
##    tt = dHKB[l]
##    rssi = tA[sc][l,:]
##    t = np.arange(len(rssi))*25.832e-3
##    plt.plot(t,rssi,label=tt)
##    #plt.title(dHKBc[sc])
##
##plt.legend()
##plt.ylim(-90,-45)
##
##
### In[56]:
##
##sc   = 3
### Sc20_S6_R2_HKBS_8
##links = [176,75,76,77,78]
##plt.figure(figsize=(20,5))
##
##for l in links:
##    tt = dHKB[l]
##    rssi = tA[sc][l,:]
##    t = np.arange(len(rssi))*25.8232e-3
##    plt.plot(t,rssi,label=tt)
##    plt.title('S20_S5_R1_HKBS 1"40')
##
##plt.legend()
##plt.ylim(-90,-45)
##
##
### In[59]:
##
##sc   = 21
### Sc20_S6_R2_HKBS_8
##links = [176,75,76,77,78]
##plt.figure(figsize=(20,5))
##
##for l in links:
##    tt = dHKB[l]
##    rssi = tA[sc][l,:]
##    t = np.arange(len(rssi))*25.832e-3
##    plt.plot(t,rssi,label=tt)
##    plt.title('S20_S5_R1_HKBS 1"40')
##
##plt.legend()
##plt.ylim(-90,-45)
##
