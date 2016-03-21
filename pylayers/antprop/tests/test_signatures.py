from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
from pylayers.antprop.signature import *
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
import pylayers.signal.bsignal as bs
from datetime import datetime
import pdb
import numpy as np
import matplotlib.pyplot as plt
import time

# Environnement
from pylayers.gis.layout import *

# Impulsion
import pylayers.signal.waveform as wvf

# Generation de la CIR
from pylayers.simul.link import *
# Sauvegarde
import cPickle as pickle
import scipy.io as sio
import time


S = Simul()
S.layout('TA-Office.ini')
AnchorNodes = {1:{'name':'Tag_1','coord':[35, 3, 1.64]},
               2:{'name':'Tag_2','coord':[39, 6, 1.64]},
               3:{'name':'Tag_2','coord':[30, 13, 1.64]},
               4:{'name':'Tag_2','coord':[22, 9, 1.64]},
               5:{'name':'Tag_2','coord':[2, 9, 1.64]},
               6:{'name':'Tag_2','coord':[2, 6, 1.64]},
               7:{'name':'Tag_2','coord':[22, 6, 1.64]},
              }


S.tx.clear()
S.rx.clear()
S.tx.filant='def.vsh3'
S.rx.filant='def.vsh3'
da ={}
dm ={}
S.tx.position
fig = plt.figure(figsize=(20,20))
fig,ax = S.L.showG('st',nodes=False, fig=fig,aw=True)
plt.axis('on')
c = (5,0)
k = AnchorNodes.keys()[c[0]]
pta = array([AnchorNodes[k]['coord'][0], AnchorNodes[k]['coord'][1], AnchorNodes[k]['coord'][2]]).reshape(3,1)
S.tx.point(pta,mode="add")
# ###### Trajectoire
S.rx.linevect(npt=60, step=0.1, ptt=[39, 10, 1.275], vec=[-1, 0, 0], mode='subst')
ps = S.rx.position[:,-1]
S.rx.linevect(npt=10, step=0.1, ptt=ps, vec=[0,-1,0], mode='append')
ps = S.rx.position[:,-1]
S.rx.linevect(npt=250, step=0.1, ptt=ps, vec=[-1,0,0], mode='append')
ps = S.rx.position[:,-1]
S.rx.linevect(npt=30, step=0.1, ptt=ps, vec=[0,-1,0], mode='append')
ps = S.rx.position[:,-1]
S.rx.linevect(npt=300, step=0.1, ptt=ps, vec=[1,0,0], mode='append')

# Waveform
wav = wvf.Waveform(typ='generic',bandGHz=2.559,fcGHz=4.9936,feGHz=100,threshdB=-10,twns=64)
fig = plt.figure(figsize=(10,10))
wav.show(fig=fig)
#fGHz = wav.bandwidth(th_ratio=100,Npt=200)
fGHz = np.linspace(2.5,7.5,60)

L=Layout('TA-Office.ini')
link = DLink(force=True,L=L, fGHz=fGHz, verbose=False)
link.fGHz=fGHz


link.a = S.tx.position[:,1]
link.b = S.rx.position[:,22]



print "eval ..."
#tic1 = time.clock()
#(ak,tauk)= link.eval(force=['sig','R','Ct','H'],alg=2015,ra_ceil_height_meter=3,ra_number_mirror_cf=1,verbose=False)
#tic2 = time.clock()
#print "Algo 2015 :",tic2-tic1,len(link.Si),len(link.Si.keys())
#(ak,tauk)= link.eval( force=['sig','R','Ct','H'],alg='20152',ra_ceil_height_meter=3,ra_number_mirror_cf=1,verbose=False)
#tic3 = time.clock()
#print "Algo 20152 :",tic3-tic2,len(link.Si),len(link.Si.keys())
#(ak,tauk)= link.eval( force=['sig','R','Ct','H'],si_algo='old',alg='5',ra_ceil_height_meter=3,ra_number_mirror_cf=1,verbose=False)
#tic4 = time.clock()
#print "Algo 5 (old propaths) :",tic4-tic3,len(link.Si),len(link.Si.keys())
#(ak,tauk)= link.eval( force=['sig','R','Ct','H'],si_algo='new',alg='5',ra_ceil_height_meter=3,ra_number_mirror_cf=1,verbose=False)
tic5 = time.clock()
#print "Algo 5 (new procone) :",tic5-tic4,len(link.Si),len(link.Si.keys())
(ak,tauk)= link.eval( force=['sig','R','Ct','H'],alg=7,ra_ceil_height_meter=3,ra_number_mirror_cf=1,verbose=False)
tic6 = time.clock()
print "Algo 7  :",tic6-tic5,len(link.Si),len(link.Si.keys())
#print "apply ..."
#ciro = link.H.applywav(wav.sfg)
#
#
#
#taumin = 160
#taumax = 165
#u = np.where((tauk>taumin) &(tauk<taumax))
#
#def display_rays(taumin,delta):
#    taumax = taumin+delta
#    u = np.where((tauk>taumin) &(tauk<taumax))
#    link.R.show(L=L,rlist=u[0],figsize=(15,8),colray='red')
#    plt.title(str(u[0]))
#    plt.figure()
#    ciroc = np.sum(ciro.y,axis=0)
#    plt.plot(ciro.x,ciroc)
#    for k in u[0]:
#        plt.plot(ciro.x,ciro.y[k,:],'b')
#    plt.xlim(taumin,taumax)
#
#
