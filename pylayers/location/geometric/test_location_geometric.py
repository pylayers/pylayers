#!/usr/bin/python
#-*- coding:Utf-8 -*-
import numpy as np

from pylayers.location.geometric.constraints.cla import *
from pylayers.location.geometric.constraints.toa import *
from pylayers.location.geometric.constraints.tdoa import *
from pylayers.location.geometric.constraints.rss import *
from pylayers.network.model import *

import matplotlib.pyplot as plt
from pylayers.location.geometric.util.cdf2 import CDF

# size of box [0, Lmax]
Lmax=20

# standard deviation
std = 5

# nb of trials
N=100

P = np.random.rand(2,N)

# anchors :
A = np.array(([[0,0],
                [0,Lmax],
                [Lmax,0],
                [Lmax,Lmax]]
                )).T

# distance / TOA
D = np.sqrt(np.sum((A[:,None,:]-P[:,:,None])**2,axis=0))
toa = D/0.3

# error 
n = std*np.random.randn(N,A.shape[1])
toa= toa + n

# TDOAs
tdoa = np.array(([toa[:,0]-toa[:,1]],[toa[:,0]-toa[:,2]],[toa[:,0]-toa[:,3]])).reshape(N,3)
# ATDOA : nb_tdoa,a1a2,xy
ATDOA= np.array(([A[:,0],A[:,1]], [A[:,0],A[:,2]], [A[:,0],A[:,3]]))




# nodes={}
# ldp={}
# ldp['TOA'] = toa
# ldp['TOA_std'] = 1.0
# ldp['TDOA'] = tdoa
# ldp['TDOA_std'] = 1.0

pe_toa = []
pe_tdoa = []
pe_toa_tdoa = []
for k in xrange(N):
    C_toa=CLA()
    C_tdoa=CLA()
    C_toa_tdoa=CLA()
    T = []
    TD = []
    # TOA
    for ia in range(A.shape[1]):
        T.append(TOA(id=0,value = toa[k,ia], std = std, p = A[:,ia]))
        C_toa.append(T[ia])
        C_toa_tdoa.append(T[ia])
    C_toa.update()
    C_toa.compute()
    pe_toa.append(C_toa.pe)

    # TDOA
    for ia in range(A.shape[1]-1):
        TD.append(TDOA(id=0,value = tdoa[k,ia], std = std, p = ATDOA[ia]) )
        C_tdoa.append(TD[ia])
        C_toa_tdoa.append(TD[ia])

    C_tdoa.update()
    C_tdoa.compute()
    pe_tdoa.append(C_toa.pe)

    C_toa_tdoa.update()
    C_toa_tdoa.compute()
    pe_toa_tdoa.append(C_toa_tdoa.pe)

pe_toa=np.array(pe_toa).T
pe_tdoa=np.array(pe_tdoa).T
pe_toa_tdoa=np.array(pe_toa_tdoa).T





RMSE={}

RMSE['toa'] = np.sqrt(np.sum((pe_toa-P)**2,axis=0))
RMSE['tdoa'] = np.sqrt(np.sum((pe_tdoa-P)**2,axis=0))
RMSE['toa+tdoa'] = np.sqrt(np.sum((pe_toa_tdoa-P)**2,axis=0))




co = ['r','g','b']
lv=[]
for ii,i in enumerate(RMSE):
    d0 = {}
    d0['values'] = RMSE[i]
    d0['bound'] = np.arange(0, 10, 0.1)
    d0['xlabel'] = 'Cumulative probablity'
    d0['ylabel'] = 'Distance [m]'
    d0['legend'] = i
    d0['markersize'] = 3
    d0['markercolor'] = 'red'
    d0['markerfrequency'] = 2
    d0['title'] = 'title'
    d0['color'] = co[ii]
    d0['marker'] = 'o'
    d0['line'] = '-'
    d0['linewidth'] = 3
    d0['filename'] = 'essai.png'
    d0['save'] = False
    lv.append(d0)
c = CDF(lv, 'fig')
c.show()


