from pylayers.measures.vna.analyzer import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
import pdb

<<<<<<< HEAD
vna = VNA()
vna.calibh5(_fileh5='mytest',
         _filecal='cal_config.ini',
=======
vna = SCPI()
_fileh5 = 'test7Dec'
# calibration
vna.calibh5(_fileh5=_fileh5,
        _filecal='cal_config.ini',
>>>>>>> af51b5f50ee897fe665933ad51b8419282bf04b8
         _filevna='vna_config.ini',
         cables=['cn1','t1','cn2'],
         author='Mamadou',
         comment='test')
<<<<<<< HEAD
vna.showcal()
=======
scanner = Scanner()
## measure
#A1 = AntArray(N=[8,1,1],dm=[0.075,0,0])
#scanner.meash5(A1, _fileh5=_fileh5,gcal=1,ical=1)
A2 = AntArray(N=[8,4,1],dm=[0.075,0.075,0])
scanner.meash5(A2, _fileh5=_fileh5,gcal=1,ical=7,Nt=4,Nr=4)
#A3 = AntArray(N=[8,8,4],dm=[0.075,0.075,0.075])
#scanner.meash5(A3, _fileh5=_fileh5,gcal=1,ical=7)
#
>>>>>>> af51b5f50ee897fe665933ad51b8419282bf04b8
