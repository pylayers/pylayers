from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *

vna = SCPI()
_filecalh5 = 'mamadou_mimo_4_8_V2'
vna.calibh5(Nr = 4,
         Nt = 8,
         _filecalh5=_filecalh5,
         _filecal='cal_config.ini',
         _filevna='vna_config.ini',
         typ='full',
         cables=['CN27','CN29','RF-OPT','OPT-RF','CN3'],
         author='Mamadou',
         comment='test MIMO calibration RF-OPT with trigger bus and with nf = 1601')