import numpy as np
import pdb
from pylayers.antprop.slab import *
from pylayers.antprop.diff import *
#
# Metalic case : MacNamara Page 202
#
x = np.logspace(-4,1,100)
y,ys,yl=FreF(x)
y2=FreF2(x)
plt.ion()
fig,ax1=plt.subplots(figsize=(10,5))
ax1.semilogx(x,np.abs(ys),'k.',label='module FreF')
ax1.semilogx(x,np.abs(yl),'k.-',label='module FreF')
ax1.semilogx(x,np.abs(y),'k',label='module FreF')
ax1.semilogx(x,np.abs(y2),'r',label='module FreF2')
plt.legend()
plt.ylabel('Magnitude',fontsize=14)
plt.ylim(0,1)
plt.xlabel('KLa',fontsize=14)
ax2=ax1.twinx()
ax2.semilogx(x,np.angle(ys)*180/np.pi,'b.',label='phase FreF')
ax2.semilogx(x,np.angle(yl)*180/np.pi,'b.-',label='phase FreF')
ax2.semilogx(x,np.angle(y)*180/np.pi,'b',label='phase FreF')
ax2.semilogx(x,np.angle(y2)*180/np.pi,'r',label='phase FreF2')
plt.ylabel('Phase degrees',fontsize=14)
plt.title(u'F(x) transition function',fontsize=14)
plt.ylim(0,50)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('F.png')
