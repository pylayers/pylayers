from pylayers.gis.layout import *
import pylayers.util.geomutil as geu
Lfile = 'scattering_nonconvex.ini'
data = '/home/niamiot/Documents/code/pylayers/data/struc/ini/'+Lfile
proj = '/home/niamiot/Documents/Pylayers_project/P1/struc/ini/'+Lfile

shutil.copyfile(data,proj)

L = Layout(Lfile,force=True)
L.build('t')
L._convex_hull()
# L.showG('st',airwalls=True)
# L._convexify()