from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
import pdb

vna = SCPI()
_filecalh5 = 'mimocal8_4_V1'
# 1. perform a new SISO calibration with the VNA
vna.calibh5(Nr = 4,
         Nt = 8,
         _filecalh5=_filecalh5,
         _filecal='cal_config.ini',
         _filevna='vna_config.ini',
         typ='full',
         cables=['CN14','CN27'],
         author='Mamadou',
         comment='first MIMO calibration RF')

# 2. Initialize the scanner
#scanner = Scanner(Nt=4,Nr=4)

#A1 = AntArray(N=[8,1,1],dm=[0.075,0,0])
#scanner.meash5(A1, _fileh5=_fileh5,gcal=1,ical=1)

# 3. Define a grid of points
#grid = AntArray(N=[10,6,2,1],dm=[0.075,0.075,0.075,0])

# 4. perform a new measurement
#
# The nature of the measurement is defined through Nt and Nr
# If Nt = 1 and Nr = 1 --> SISO case
# Otherwise -> MISO , SIMO or MIMO case.
#

#scanner.meash5(grid,_fileh5=_fileh5,gcal=1,ical=7)
#A3 = AntArray(N=[8,8,4],dm=[0.075,0.075,0.075])
#scanner.meash5(A3, _fileh5=_fileh5,gcal=1,ical=7)
#
