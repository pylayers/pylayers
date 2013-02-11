#!/usr/bin/python
#-*- coding:Utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from pylayers.gis.layout import *
from pylayers.antprop.signature import *

# load the layout graphs
def showr2(L,r2d,tx,rx,k,l):
    col = ['r','b','g','c','m','k','y']
    r = r2d[str(k)]
    pts = r['pt']
    sig = r['sig']
    fig,ax = showsig(L,sig[:,:,l],tx,rx) 
    sh = np.shape(pts)
    x = np.hstack((tx[0],pts[0,:,l],rx[0]))
    y = np.hstack((tx[1],pts[1,:,l],rx[1]))
    plt.plot(x,y,col[k])
    plt.title(sig[:,:,l])
    return fig,ax

def showr2d(L,r2d,tx,rx):
    """
    r2d['pt'] : nd,ni,nr
    """
    L.display['thin']=True
    col = ['r','b','g','c','m','k','y']
    fig,ax = L.showGs()
    for k in r2d:
        r = r2d[k]
        pts = r['pt']
        sh = np.shape(pts)
        for r in range(sh[2]):
            x = np.hstack((tx[0],pts[0,:,r],rx[0]))
            y = np.hstack((tx[1],pts[1,:,r],rx[1]))
            plt.plot(x,y,col[eval(k)])
    return fig,ax

def showsig(L,s,tx,rx):
    L.display['thin']=True
    fig,ax = L.showGs()
    L.display['thin']=False
    L.display['edlabel']=True
    fig,ax = L.showGs(fig=fig,ax=ax,edlist=s[0,:],width=4)
    plt.plot(tx[0],tx[1],'x')
    plt.plot(rx[0],rx[1],'+')
    plt.title(str(s[0,:])+str(s[1,:]))
    L.display['edlabel']=False
    return fig,ax

#strucname = 'TA-Office'
strucname = 'DLR'
L = Layout(strucname+'.ini')
L.boundary()
print L.ax
#try:
#    L.dumpr()
#except:
L.build()
L.dumpw()
#tx = np.array([8., 8., 1.])
#rx = np.array([30., 11., 2.])
tx = np.array([1., 0., 1.])
rx = np.array([8., 0., 2.])

#L = Layout('TA-Office.str')
#L.build()
#tx = np.array([20, 8, 1])
#rx = np.array([35, 6, 2])


S = Signatures(L, tx, rx)

print "Calcul signatures"
#s1 = S.get_sigslist(tx, rx)
s1 = S.run(tx,rx,5)
print "Fin calcul signatures"

#print "signatures --> rayons "
#r2d = S.sigs2rays(s1)
r2d = S.rays(s1)
##print "fin signatures --> rayons "
##
#r22 = r2d['2']
#pt2 = r22['pt']
#sig2 = r22['sig']
#pt2 = np.swapaxes(pt2,0,2)
#pt2 = np.swapaxes(pt2,1,2)
#tx2 = np.kron(np.ones(2),tx).reshape(2,3,1)
#rx2 = np.kron(np.ones(2),rx).reshape(2,3,1)
#tx2[:,2,:]=0
#rx2[:,2,:]=0
#pt  = np.concatenate((tx2,pt2,rx2),axis=2)
#vsi = pt[:, :, 1:] - pt[:,:,:-1]
#si  = np.sqrt(np.sum(vsi*vsi, axis=1))
#alpha = np.cumsum(si,axis=1)
#c  = alpha[:,-1].reshape(2,1)
#alpha = alpha/c
#pt[:,2,1:]= tx[2]+alpha*(rx[2]-tx[2])
#
#
showr2d(L,r2d,tx,rx)
print "rayons 2D --> rayons3D "
#rays3d = S.ray2D3D(r2d)
#print "fin rayons 2D --> rayons3D "
##
#S.show3(rays=rays3d,strucname=strucname)
##
##
##
#s    = np.array([[5,1,8],[1,1,2]])
#sig  = Signature(s)
#rsig = sig.sig2ray(L,tx[0:2],rx[0:2])
#sig.ev(L)
#M = sig.image(tx[0:2])
#Y = sig.backtrace(tx[0:2],rx[0:2],M)
#plt.plot(M[0,:],M[1,:],'ob')
#plt.plot(Y[0,:],Y[1,:],'xk')
#fig,ax = showr2(L,r2d,tx[0:2],rx[0:2],3,4)
#plt.show()
#room8 = L.Gt.node[8]
#polyg8 = room8['polyg']
#vnodes8 = room8['vnodes']
#udeg1 = []
#udeg2 = []
#for ik, inode in enumerate(vnodes8):
#    deg = L.Gs.degree(inode)
#    if vnodes8[0] < 0:
#        index = ik / 2
#    else:
#        index = (ik - 1) / 2
#    if inode < 0:
#        if deg == 2:
#            udeg2.append(index)
#        if deg == 1:
#            udeg1.append(index)    # warning not used
#Gv = polyg8.buildGv(show=True,udeg2=udeg2)
#L.showGs()
#nx.draw_networkx_edges(L.dGv[8],L.Gs.pos,nx.edges(L.dGv[8],nbunch=[47]))
