import numpy as np
import matplotlib.pyplot as plt
from pylayers.antprop.slab import *
#
# get the slab database
#
sl  = SlabDB('matDB.ini','slabDB.ini')
#
## define a material and a thickness
##
#lname  = ['BRICK']
#lthick = [0.5]
#sl.add('test',lname,lthick)
##
## define a time base
##
df = 0.01
fGHz   = np.arange(1,5,df)
##
## define an angular value
##
theta = 1
#plt.ion()
#sl['test'].ev(fGHz,theta,compensate=True)
#sl['test'].plotwrtf(typ='mod')
#plt.figure()
#sl['test'].plotwrtf(typ='angle')
#plt.show()
#T   = sl['test'].T
#To  = T[:,0,0,0]
#Tp  = T[:,0,1,1]
#ao  = np.unwrap(np.angle(To))
#ap  = np.unwrap(np.angle(Tp))
#delayo = np.diff(ao)/(2*np.pi*df)
#delayp = np.diff(ap)/(2*np.pi*df)
#plt.figure()
#plt.plot(fGHz[0:-1],-delayo*0.3*100)
#plt.plot(fGHz[0:-1],-delayp*0.3*100)
#plt.xlabel('frequency(GHz)')
#plt.ylabel('excess distance (cm)')
##    #    thick  = [0.50,0]
##    thick2  = [0.1,0.3,0.1]
##    theta  = np.arange(0,np.pi/2,0.01,dtype=np.float64)
##    fGHz   = np.arange(0.4,2.3,0.05)
##    #    theta  = array([np.pi/2-0.01])
##    #    theta  = array([np.pi/4])
##    #    fGHz   = array([5.0])
##    #    fGHz      = array([2.4])
##    #    fGHz      = array([0.3])
##    S1     = MLayer(lmat,thick2,fGHz,theta)
##    S1.RT()
##
##    plt.figure()
##    S1.plotwrta(0)
##    plt.figure()
##    S1.plotwrtf(0)
##    #    S1 = MLayer(lmat,thick,fGHz,theta)
##    S1.pcolor()
##    #    S1.plotwrta()
##    #    S2 = MLayer(lmat2,thick2,fGHz,theta)
#
##    #    S2.plotwrta()
##    #    II.plotwrta(0)
##    #    II.plotwrtf(0)
##    #    II.pcolor(f,theta)
##    #
##
##    #    sl['WALL'].eval(f,theta)
#sl['WOOD'].ev(fGHz,theta,compensate=True)
#sl['WOOD'].plotwrtf(typ='mod')
sl['METAL'].ev(fGHz,theta,compensate=False)
sl['METAL'].plotwrtf(typ='mod')
#    #
#    #
#    #
#    #       L=sl['WALL'].loss0(f)
#    #    Lo,Lp=sl['WALL'].losst(f,pi/3)
