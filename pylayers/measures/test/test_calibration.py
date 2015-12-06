from pylayers.measures.vna.analyzer import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *

vna = VNA()
vna.calibh5(_fileh5='mytest',
         _filecal='cal_config.ini',
         _filevna='vna_config.ini',
         cables=['cn1','t1','cn2'],
         author='Mamadou',
         comment='test')
vna.showcal()
