from pylayers.measures.vna.E5072A import *
from pylayers.measures.parker.smparker import *
from pylayers.antprop.channel import *
from pylayers.antprop.aarray import *
import pylayers.util.pyutil as pyu

#Nr = input ("Number of receiving antennas : ")
#Nt = input ("Number of transmitting antennas : ")
Nr = 2
Nt = 2
author = 'bu'
comment = 'no comment'
# Full MIMO calibration file (without extension .h5)
_filecalh5 = 'calMIMO'
# MIMO measurement file
_filemesh5 = 'mesMIMO'

emulated = True
mimocal = False

filecalh5 = pyu.getlong(_filecalh5,pstruc['DIRMES'])+'.h5'
if os.path.isfile(filecalh5):
    C = Mesh5(filecalh5)
    if len(C.gcal.keys())>1:
        print " The Full MIMO calibration file contains the following group of calibration (gcal) : "
        for g in C.gcal.keys():
            print g.replace('cal','')
        num = raw_input("Select group of mimo calibration number ")
        gcalm = num   # group calibration mimo
        gcal = 'cal'+ gcalm
    else:
        gcalm = '1'
        gcal = 'cal1'
    print "The Full MIMO calibration group :", gcal, " , contains the following configurations"
    print "-----------------------------------------------------------------------\n"
    for k1 in C.gcal[gcal]:
        print "ical : ",k1
        for k2 in C.gcal[gcal][k1]:
            if k2 in ['time','ifbhz','nf']:
                print " ",k2,':',C.gcal[gcal][k1][k2]
else:
    print "Warning there is no full MIMO calibration file in your measurement directory"
    rep  = raw_input("Do you want to create one  ? (Y/N) ")
    if rep=='Y':
        mimocal=True
        gcalm = '1'
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
         _filecal='calibration_config.ini',
         _filevna='vna_configuration.ini',
         typ='full',
         cables = ['CN27','CN29','RF-OPT','OPT-RF','CN3'],
         author = author,
         comment = comment)
   vna.close()

# check if measurement file exists
filemesh5 = pyu.getlong(_filemesh5,pstruc['DIRMES'])+'.h5'

if os.path.isfile(filemesh5):
    M = Mesh5(filemesh5)
    lkeys = M.gcal.keys()
    lcal = [ c for c in lkeys if 'cal' in c ]
    
    print "\n\nThe measurement file contains the following group of calibration (gcal) : "
    for g in lcal:
        print g.replace('cal','')
    num = raw_input("Select group of single calibration number or (N for a new one) :")
    if (num == 'N') or (num =='n'):
        calibration = True
        measure = False
    else: 
        calibration = False
        measure = True
        gcals = int(num)    # group of calibration single 
        gcal = 'cal'+ str(gcals)

    print "The calibration group :", gcal, " , contains the following configurations"
    print "-----------------------------------------------------------------------\n"
    for k1 in M.gcal[gcal]:
        print "ical : ",k1
        for k2 in M.gcal[gcal][k1]:
            if k2 in ['time','ifbhz','nf','power']:
                print " ",k2,':',M.gcal[gcal][k1][k2]
    sical = raw_input("Select a configuration number for the next measurement series:")
    ical = int(sical)

else:
    print "Warning there is no measurement file in your measurement directory"
    rep  = raw_input("Do you want to make a single calibration (Y/N) ")
    if rep=='Y':
        calibration=True
    else:
        calibration=False
        exit
    
if calibration:
   print "Phase 2 : MIMO single channel calibration"
   #author = input('author : ')
   #comment  = input('comment : ')
   vna = SCPI(emulated=emulated)
   vna.calibh5(Nr = Nr,
            Nt = Nt,
            _filecalh5=_filecalh5,
            _filemesh5=_filemesh5,
            _filecal='calibration_config.ini',
            _filevna='vna_configuration.ini',
            typ='single',
            gcalm = gcalm,
            cables=['CN27','CN29','RF-OPT','OPT-RF','CN3'],
            author=author,
            comment=comment)
   vna.close()
if measure:
   print "Phase 3 : MIMO single channel calibration"
   print  "Measure started !"
   sheight = raw_input('Transmitting antenna height (meters) : ')
   height = float(sheight)
   comment  = input('comment : ')
   # 2. Initialize the scanner
   scanner = Scanner(Nr=Nr,Nt=Nt,emulated=emulated)
   #A1 = AntArray(N=[20,1,1,1],dm=[0.01,0,0])

   A1 = AntArray(N=[5,1,1,1],max=[-0.35,0,0,0],min=[0.35,0,0,0],mode='grid')
   tic = time.time()
   scanner.meas(A1,
              _fileh5=_filemesh5,
              gcal = gcals,
              ical = ical,
              vel = 15,
              Nmeas = 1,
              pAnt = np.array([1.6,5.2,height]),
              vAnt = np.array([1.0,0.0,0.0]),
              comment = comment,
              author = 'mamadou',
              )
   toc = time.time()
   print "End of measurement time (s) :",toc-tic
