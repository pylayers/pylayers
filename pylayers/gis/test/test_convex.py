from pylayers.gis.layout import *
import pylayers.util.geomutil as geu
import os

Lfile = 'WHERE1_nc.ini'
Lfile = 'TA-Office.ini'


data = os.environ['PYLAYERS']+'/data/struc/ini/'+Lfile
proj = os.environ['BASENAME']+'/struc/ini/'+Lfile

shutil.copyfile(data,proj)

L = Layout(Lfile,force=True)
L.build('t',convex=True)
# L._convex_hull()
# # L.showG('st',airwalls=True)
# L._convexify()


# p1 = geu.Polygon(L.ax,delta=5)
# L.ma = L.mask()
# p2 = p1.difference(L.ma)
# boundary = geu.Polygon(p2)
# boundary.vnodes = L.ma.vnodes
# L.Gt.add_node(0,polyg=boundary)
# L.Gt.add_node(0, indoor = False)
# L.Gt.pos[0]=(L.ax[0],L.ax[2])


# # all segments of the building boundary
# nseg = filter(lambda x : x >0 , boundary.vnodes)
# # air segments of the building boundary
# nsegair = filter(lambda x : x in L.name['AIR'],nseg)
# # wall segments of the building boundary
# nsegwall = filter(lambda x : x not in L.name['AIR'],nseg)

# #
# # ldiffin  : list of indoor diffraction points
# # ldiffout : list of outdoor diffraction points (belong to layout boundary)
# #

# L.ldiffin  = filter(lambda x : x not in boundary.vnodes,L.ldiff)
# L.ldiffout = filter(lambda x : x in boundary.vnodes,L.ldiff)
