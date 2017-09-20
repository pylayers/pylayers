# -*- coding: utf-8 -*-
r"""
=========================================
Indoor Radio Coverage with Motley Keenan
=========================================

This is an example of Motley-Keenan coverage 

The default configuration file coverage.ini. 
The antenna assoaciated to each Access Point are handled in the module 
pylayers/antprop/antenna.py

    [grid]
    nx = 80
    ny = 40
    boundary = [20,0,30,20]
    mode = full ; file, zone , full 
    file = 'points.ini'

    [layout]
    filename = TA-Office.lay ; 0 40 0 15
    ;filename = W2PTIN.lay
    ;filename = Lstruc.lay

    [ap]
    0 = {'name':'room1','wstd':'ieee80211b','p':(1,12,1.2),'PtdBm':0,'chan':[11],'on':True,'ant':'Gauss','phideg':90} 
    1 = {'name':'room2','wstd':'ieee80211b','p':(10,2,1.2),'PtdBm':0,'chan':[11],'on':True,'ant':'Omni','phideg':0} 
    2 = {'name':'room3','wstd':'ieee80211b','p':(20,1,1.2),'PtdBm':0,'chan':[11],'on':True,'ant':'Gauss','phideg':75} 
    3 = {'name':'room4','wstd':'ieee80211b','p':(36.5,1.5,1.2),'PtdBm':0,'chan':[11],'on':True,'ant':'Gauss','phideg':120} 
    4 = {'name':'room5','wstd':'ieee80211b','p':(25,12,1.2),'PtdBm':0,'chan':[11],'on':True,'ant':'Gauss','phideg':180} 

    [rx]
    temperaturek = 300
    noisefactordb = 0 

    [show]
    show = True

"""
from pylayers.antprop.coverage import *
import matplotlib.pyplot as plt
import time
import pdb
C = Coverage()
C.L._filename
C.tx = np.array((39,1))
start = time.time()
C.cover()
finish = time.time()
C.show(figsize=(16,9))
#print 'Tx position: ',C.tx 
#print 'Coverage in %1.2f seconds' % (finish-start)
plt.show()
