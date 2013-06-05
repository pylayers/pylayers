from pylayers.gis.furniture import *
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
import matplotlib.pylab as plt
import numpy as np

_filefur = 'Furw1.ini'
filefur = pyu.getlong(_filefur, pstruc['DIRSTRUC'])
config = ConfigParser.ConfigParser()
config.read(filefur)
furname = config.sections()
lfur = []
fig =plt.figure()
ax =fig.gca()

pt = np.array([-20,10])
ph = np.array([-10,14])
plt.plot([pt[0],ph[0]],[pt[1],ph[1]])
lpt =[]
for name in furname:
    F = Furniture()
    F.load('Furw1.ini',name)
    if F.Matname == 'METAL':
        lfur.append(F)
        pos = F.position()
        lpt.append(pos)
        p   = np.array(pos).T
        d1,d2,h = geu.dptseg(p,pt,ph)
        if ((np.prod(d1>0)) & (np.prod(d2>0))):
            print d1
            print d2 
            print h 
            fig,ax=F.show(fig,ax)
apt =np.array(lpt)
sh =np.shape(apt)
nfur = sh[0]

