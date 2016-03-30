from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
from pylayers.antprop.aarray import *
import pylayers.util.pyutil as pyu

#Nr = input ("Number of receiving antennas : ")
#Nt = input ("Number of transmitting antennas : ")
Nr = 4
Nt = 8
author = 'bu'
comment = 'no comment'
# Full MIMO calibration file
_filecalh5 = 'calMIMO.h5'
# MIMO measurement file
_filemesh5 = 'mesMIMO.h5'

mimocal = False
emulated = True
calibration = False
measure = False

filecalh5 = pyu.getlong(_filecalh5,pstruc['DIRMES'])
if os.path.isfile(filecalh5):
    C = Mesh5(filecalh5)
    if len(C.gcal.keys())>1:
    	print " The Full MIMO calibration file contains the following group of calibration (gcal) : "
    	for g in C.gcal.keys():
    		print g.replace('cal','')
    	num = raw_input("Select group of calibration number ")
    	gcal = 'cal'+int(num)
    else:
    	gcal = 'cal1'
    print "The calibration group :", gcal, " , contains the following configurations"
    print "-----------------------------------------------------------------------\n"
    for k1 in C.gcal[gcal]:
    	print k1
    	for k2 in C.gcal[gcal][k1]:
    		if k2 in ['time','ifbhz','nf']:
    			print " ",k2,':',C.gcal[gcal][k1][k2]
else:
    print "Warning there is no full MIMO calibration file in your measurement directory"
    rep  = raw_input("Do you want to create one  ? (Y/N) ")
    if rep=='Y':
        mimocal=True
    else:
        mimocal=False




if mimocal:
   print "Phase 1: MIMO full calibration"
   #author = input('author : ')
   #comment  = input('comment : ')
   vna = SCPI(emulated=emulated)
   vna.calibh5(Nr = Nr,
         Nt = Nt,
         _filecalh5=_filecalh5,
         _filecal='cal_configuration.ini',
         _filevna='vna_configuration.ini',
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
