from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
from pylayers.antprop.aarray import *


#
# _filecalh5 = 'MIMO_8_4_14janv_V1'
# vna.calibh5(Nr = 4,
#          Nt = 8,
#          _filecalh5=_filecalh5,
#          _filecal='cal_config.ini',
#          _filevna='vna_config.ini',
#          typ='full',
#          cables=['CN14','CN27'],
#          author='Mamadou',
#          comment='test MIMO calibration RF without nf=1601')


Nr = input ("Number of receiving antennas : ")
Nt = input ("Number of transmitting antennas : ")
author = 'bu'
comment = 'no comment'
_filecalh5 = 'calMIMO_'+str(Nr)+'_'+str(Nt)
_filemesh5 = 'mesMIMO_'+str(Nr)+'_'+str(Nt)

emulated = True
mimocal = True
calibration = True
measure = True

if mimocal:
    print "Phase 1: MIMO full calibration"
    #author = input('author : ')
    #comment  = input('comment : ')
    vna = SCPI(emulated=emulated)
    vna.calibh5(Nr = Nr,
          Nt = Nt,
          _filecalh5=_filecalh5,
          _filecal='cal_config.ini',
          _filevna='vna_config.ini',
          typ='full',
          cables = ['CN27','CN29','RF-OPT','OPT-RF','CN3'],
          author = author,
          comment = comment)
    vna.close()


if calibration:
    print "Phase 2 : MIMO single channel calibration"
    #author = input('author : ')
    #comment  = input('comment : ')
    vna = SCPI(emulated=emulated)
    vna.calibh5(Nr = Nr,
             Nt = Nt,
             _filecalh5=_filecalh5,
             _filemesh5=_filemesh5,
             _filecal='cal_config.ini',
             _filevna='vna_config.ini',
             typ='single',
             gcalm = '1',
             cables=['CN27','CN29','RF-OPT','OPT-RF','CN3'],
             author=author,
             comment=comment)
    vna.close()
if measure:
    print "Phase 3 : MIMO single channel calibration"
    print  "Measure started !"
    #author = input('author : ')
    #comment  = input('comment : ')
    # 2. Initialize the scanner
    scanner = Scanner(Nr=Nr,Nt=Nt)
    #A1 = AntArray(N=[20,1,1,1],dm=[0.01,0,0])

    A1 = AntArray(N=[1000,1,1,1],max=[0.35,0,0,0],min=[-0.35,0,0,0],mode='grid')
    tic = time.time()
    scanner.meas(A1,
               _fileh5=_filemesh5,
               gcal = 1,
               ical = 3,
               vel = 15,
               Nmeas = 1,
               pAnt = np.array([1.6,5.2,1.6]),
               vAnt = np.array([1.0,0.0,0.0]),
               comment = 'first measurement of MIMO 8X4 70cm over axis X',
               author = 'mamadou',
               )
    toc = time.time()
    print "End of measurement time (s) :",toc-tic
