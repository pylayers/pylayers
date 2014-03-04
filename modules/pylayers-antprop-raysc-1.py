from pylayers.antprop.raysc import *
from pylayers.util.project import *
from pylayers.gis.layout import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
g = GrRay3D()
lfile = g.dir()
n = len(lfile)
file0  = lfile[0]
s1 = file0.split('_')
_filestr = s1[0]+'.str'
L = Layout()
L.load(_filestr)
f,a = L.showGs()
g.load(file0,L)
f1 = g.show(f,a,np.arange(10))
plt.show()
f,a = L.showGs()
f2 = g.show(f,a,np.arange(100))
plt.show()
f,a = L.showGs()
f3 = g.show(f,a,np.arange(300))
plt.show()
