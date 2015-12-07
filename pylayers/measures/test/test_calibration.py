from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
import pdb

vna = SCPI()
_fileh5 = 'test7Dec'
# calibration
vna.calibh5(_fileh5=_fileh5,
        _filecal='cal_config.ini',
         _filevna='vna_config.ini',
         cables=['CN1','T','CN3'],
         author='Mamadou',
         comment='test')
scanner = Scanner()
## measure
#A1 = AntArray(N=[8,1,1],dm=[0.075,0,0])
#scanner.meash5(A1, _fileh5=_fileh5,gcal=1,ical=1)
A2 = AntArray(N=[8,4,1],dm=[0.075,0.075,0])
scanner.meash5(A2, _fileh5=_fileh5,gcal=1,ical=7)
#A3 = AntArray(N=[8,8,4],dm=[0.075,0.075,0.075])
#scanner.meash5(A3, _fileh5=_fileh5,gcal=1,ical=7)
#
