from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *

vna = SCPI()
_fileh5 = 'test4Dec'
# calibration
vna.calibh5(_fileh5='test4Dec',
        _filecal='cal_config.ini',
         _filevna='vna_config.ini',
         cables=['CN1','T','CN3'],
         author='Mamadou',
         comment='test')
scanner = Scanner()
# measure
A = AntArray()
scanner.meash5(A, _fileh5=_fileh5,gcal=1,ical=1)

