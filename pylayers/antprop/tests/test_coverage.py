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
C2 = Coverage('coverage2.ini')
C2.cover()
C2.show(typ='loss',vmin=-90,vmax=-20)


