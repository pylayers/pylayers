import pylayers.gis.ezone as ez
from pylayers.gis.gisutil import ent,ext2qt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import smopy
from cartopy import config
import cartopy.crs as ccrs
fig = plt.figure(figsize=(12,12))
white = np.zeros((10,10))
ax = fig.add_subplot(111)
z = ez.Ezone('N48W002')
z.loadh5()
z.rebase()
zoom=11
p = (48.721095,-1.830548)
print "p : ",p
xtile,ytile=smopy.deg2num(p[0],p[1],zoom,do_round=True)
print "xtile,ytile : ",xtile,ytile
(lat0,lon0)=smopy.num2deg(xtile,ytile,zoom,do_round=True)
(lat1,lon1)=smopy.num2deg(xtile+1,ytile+1,zoom,do_round=True)
print "lat,lon WN",lat0,lon0
print "lat,lon ES",lat1,lon1

#mp = smopy.Map((lat1,lon0,lat0,lon1),z=zoom)
mp = smopy.Map((48,-2,49,-1),z=zoom)
##f,a = z.show(alpha=0.3)
box_tile = mp.box_tile
print box_tile
L_ll,l_ll=smopy.num2deg(box_tile[0],box_tile[1]+1,zoom)
L_ur,l_ur=smopy.num2deg(box_tile[2]+1,box_tile[3],zoom)
extent_true = np.array((l_ll,l_ur,L_ll,L_ur))
print extent_true
#print extent_true
##print z.extent
f,a = z.show(fig=fig,ax=ax,alpha=0.4)
#f,a=plt.subplots(1,1)
im1 = a.imshow(mp.img,extent=extent_true,alpha=0.6)
im2 = a.imshow(white,extent=(-2.2,-0.9,47.9,49.1),alpha=0)
a.plot(p[1],p[0],'ob')
###mp.box_tile=(0,0,73000,111000)
###mp.h=73000
###mp.w=111000
###mp.box_tile=(0,111000,73000,0)
###mp.xmin = 0
###mp.ymin=0
###ax = mp.show_mpl(figsize=(20,10),alpha=1)
##fig=plt.gcf()
###z.extent_c=(0,1024,0,1280)
###z.extent_c=(506,509,351,355)
###print z.extent_c
a = z.cover(Ht=2,Hr=2,Rmax=10000)
##
