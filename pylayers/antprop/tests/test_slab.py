from pylayers.antprop.slab import *
sl  = SlabDB('matDB.ini','slabDB.ini')
#    #sl.choose()
#    print "sl  : SlabDB object"
##    mat.save2(fileout)
#    mat = sl.mat
#    m1=mat['AIR']
#    #m2=mat['METAL']
#    #m2=mat['PLATRE-57GHz']
#    m2=mat['PLASTER']
#    #    m2=mat['CONCRETE']
#    #    m2.epr=6.25
#    #    m2.sigma=0
#    #    m2=Mat(name='TEST',epr=2.25)
#
#    #lmat   = [m1,m2,m1,m2,m1]
#    lmat   = [m2,m1,m2]
#    #    lmat   = [m1,m2,m1]
#    #    lmat   = [m2]
#    #    thick  = [0.15]
#    #    thick  = [0.50,0]
#    thick2  = [0.1,0.3,0.1]
#    theta  = np.arange(0,np.pi/2,0.01,dtype=np.float64)
#    fGHz   = np.arange(0.4,2.3,0.05)
#    #    theta  = array([np.pi/2-0.01])
#    #    theta  = array([np.pi/4])
#    #    fGHz   = array([5.0])
#    #    fGHz      = array([2.4])
#    #    fGHz      = array([0.3])
#    S1     = MLayer(lmat,thick2,fGHz,theta)
#    S1.RT()
#
#    plt.figure()
#    S1.plotwrta(0)
#    plt.figure()
#    S1.plotwrtf(0)
#    #    S1 = MLayer(lmat,thick,fGHz,theta)
#    S1.pcolor()
#    #    S1.plotwrta()
#    #    S2 = MLayer(lmat2,thick2,fGHz,theta)

#    #    S2.plotwrta()
#    #    II.plotwrta(0)
#    #    II.plotwrtf(0)
#    #    II.pcolor(f,theta)
#    #
#
#    #    sl['WALL'].eval(f,theta)
#    #    sl['METAL'].eval(1.0,0)
#    #
#    #
#    #
#    #       L=sl['WALL'].loss0(f)
#    #    Lo,Lp=sl['WALL'].losst(f,pi/3)
