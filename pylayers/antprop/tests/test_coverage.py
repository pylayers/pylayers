import time
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.util.utilnet import str2bool
from pylayers.gis.layout import Layout
from pylayers.antprop.coverage import *
from pylayers.network.model import *


C1 = Coverage('coverage.ini')
C1.cover()
C1.show(typ='loss',vmin=-80,vmax=-20)
plt.savefig('test-coverage-1.png')
C2 = Coverage('coverage2.ini')
C2.cover()
C2.show(typ='loss',vmin=-90,vmax=-20)
plt.savefig('test-coverage-2.png')

C=Coverage('where1cov.ini')
C.cover()
fig =plt.figure(figsize=(18,6))
ax=fig.add_subplot(221)
f,a=C.show(fig=fig,ax=ax,typ='pr',a=0)
ax=fig.add_subplot(222)
f,a=C.show(fig=fig,ax=ax,typ='pr',a=1)
ax=fig.add_subplot(223)
f,a=C.show(fig=fig,ax=ax,typ='pr',a=2)
ax=fig.add_subplot(224)
f,a=C.show(fig=fig,ax=ax,typ='pr',a=3)
plt.savefig('test-coverage-3.png')

