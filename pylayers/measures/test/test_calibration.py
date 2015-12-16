from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
import pdb

vna = SCPI()
_fileh5 = 'test15Dec'
# 1. perform a new SISO calibration with the VNA
vna.calibh5(_fileh5=_fileh5,
        _filecal='cal_config.ini',
         _filevna='vna_config.ini',
         cables=['CN1','T','CN3'],
         author='Mamadou',
         comment='test')

# 2. Initialize the scanner
scanner = Scanner(Nt=1,Nr=1)

#A1 = AntArray(N=[8,1,1],dm=[0.075,0,0])
#scanner.meash5(A1, _fileh5=_fileh5,gcal=1,ical=1)

# 3. Define a grid of points
grid = AntArray(N=[10,6,2,1],dm=[0.075,0.075,0.075,0])

# 4. perform a new measurement
#
# The nature of the measurement is defined through Nt and Nr
# If Nt = 1 and Nr = 1 --> SISO case
# Otherwise -> MISO , SIMO or MIMO case.
#

scanner.meash5(grid,_fileh5=_fileh5,gcal=1,ical=7)
#A3 = AntArray(N=[8,8,4],dm=[0.075,0.075,0.075])
#scanner.meash5(A3, _fileh5=_fileh5,gcal=1,ical=7)
#