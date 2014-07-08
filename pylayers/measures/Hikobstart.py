import os
import pdb
import sys
import pandas as pd
import numpy as np
import scipy.io as io
from pylayers.util.project import *
import matplotlib.pyplot as plt


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

CORMORAN='/media/niamiot/DONNEES/svn2/measures/CORMORAN'
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
    for k in range(1,17):
        datao[k]={}
        for l in range(1,17):
            if k!=l:
                dt = data2[k][l]
                datao[k][l] = dt[tstart[k-1][l-1]:]
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


# print "get serie 12"
# # TO BE DONE
# tstart12 = matstart(data12,1).astype('int')
# tstart12[1,6] = 4332
# tstart12[8,6] = 4209
# tstart12[10,6] = 4289
# s,data13 = extracts(data12,tstart12)

# S12 = pd.DataFrame(s)


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
