from pylayers.coverage import *
import time
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.util.utilnet import str2bool
from pylayers.gis.layout import Layout
from pylayers.antprop.multiwall import *
from pylayers.antprop.coverage import *
from pylayers.network.model import *


L = Layout('TA-Office.ini')
L.dumpr()
A=np.array((4,1)) # defining transmitter position 
B=np.array((30,12)) # defining receiver position
fGHz = 2.4
r = np.array((B,B))

