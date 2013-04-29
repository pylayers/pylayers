from pylayers.gis.layout import *
from pylayers.antprop.signature import *
from pylayers.antprop.channel import *
import pylayers.signal.waveform as wvf
import networkx as nx
import numpy as np
import time
import logging

L = Layout('WHERE1_clean.ini')
#L = Layout('defstr2.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
#L.build()
#L.dumpw()
#L.buildGi()
nc1 = 6#5
nc2 = 25#37

poly1 = L.Gt.node[nc1]['polyg']
cp1 = poly1.centroid.xy

poly2 = L.Gt.node[nc2]['polyg']
cp2 = poly2.centroid.xy
ptx = np.array([cp1[0][0],cp1[1][0],1.5])
prx = np.array([cp2[0][0]+0.5,cp2[1][0]+0.5,1.5])
print ptx
print prx
d = np.sqrt(np.dot((ptx-prx),(ptx-prx)))
tau = d/0.3
print d,tau



logging.info('Signature')
S = Signatures(L,nc1,nc2)
a =time.time()
logging.info('Calculate signature')
#S.run2(cutoff=6,dcut=3)
S.run(cutoff=2)
b=time.time()
print b-a

for i in L.Gi.nodes():
    ei = eval(i)
    if type(ei)!= int:
        if ei[0] == 354:
            print i

#Gsi.add_node('Tx')
#Gsi.pos['Tx']=tuple(ptx[:2])

#for i in L.Gt.node[nc1]['inter']:
#    if i in  Gsi.nodes():
#        Gsi.add_edge('Tx',i)

#Gsi.add_node('Rx')
#Gsi.pos['Rx']=tuple(prx[:2])

#for i in L.Gt.node[nc2]['inter']:
#    if i in  Gsi.nodes():
#        Gsi.add_edge(i,'Rx')

#print 'signatures'
#co = nx.dijkstra_path_length(Gsi,'Tx','Rx')
#sig=list(nx.all_simple_paths(Gsi,'Tx','Rx',cutoff=co+2))



#b=time.time()
#print b-a
#f,ax=L.showG('t')
#nx.draw(Gsi,Gsi.pos,ax=ax)
#plt.show()

##S.run(L,metasig,cutoff=3)
#print "r = S.rays "
r = S.rays(ptx,prx)
print "r3 = r.to3D "
r3 = r.to3D()
print "r3.locbas "
r3.locbas(L)
#print "r3.fillinter "
r3.fillinter(L)
r3.show(L)
plt.show()
##
#config = ConfigParser.ConfigParser()
#_filesimul = 'default.ini'
#filesimul = pyu.getlong(_filesimul, "ini")
#config.read(filesimul)
#fGHz = np.linspace(eval(config.get("frequency", "fghzmin")),
#                   eval(config.get("frequency", "fghzmax")),
#                   eval(config.get("frequency", "nf")))
#
#Cn=r3.eval(fGHz)
#
#Cn.freq=Cn.fGHz
#sco=Cn.prop2tran(a='theta',b='theta')
#wav = wvf.Waveform()
#ciro = sco.applywavB(wav.sfg)
#
##raynumber = 4
#
##fig=plt.figure('Cpp')
##f,ax=Cn.Cpp.plot(fig=fig,iy=np.array(([raynumber])))
#
##r3d.info(raynumber)
## plt.show()
##
##
##
###
###c11 = r3d.Ctilde[:,:,0,0]
###c12 = r3d.Ctilde[:,:,0,1]
###c21 = r3d.Ctilde[:,:,1,0]
###c22 = r3d.Ctilde[:,:,1,1]
###
###
###
###Cn=Ctilde()
###Cn.Cpp = bs.FUsignal(r3d.I.f, c11)
###Cn.Ctp = bs.FUsignal(r3d.I.f, c12)
###Cn.Cpt = bs.FUsignal(r3d.I.f, c21)
###Cn.Ctt = bs.FUsignal(r3d.I.f, c22)
###Cn.nfreq = r3d.I.nf
###Cn.nray = r3d.nray
###Cn.tauk=r3d.delays
###
###raynumber = 4
###
###fig=plt.figure('Cpp')
###f,ax=Cn.Cpp.plot(fig=fig,iy=np.array(([raynumber])))
###
##
##
##
##
##
##
