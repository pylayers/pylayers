from pylayers.simul.simulem import *
import ConfigParser
import pylayers.util.pyutil as pyu
import itertools


# Load simulnet_data configuration file
simcfg = ConfigParser.ConfigParser()
simcfg.read(pyu.getlong('pyray.ini',pstruc['DIRNETSAVE']))

Lfilename = simcfg.get('layout','layoutname')
AG = eval(simcfg.get('nodes','ag'))
AP = eval(simcfg.get('nodes','ap'))
uptime = eval(simcfg.get('simulation','updatetime'))


# create a Simul object with the correct layout

S = Simul()
S.layout(Lfilename,'matDB.ini','slabDB.ini')    
S.clean_project(verbose=True)

#
### STEP 1 : all mobile node with all agent
#

lap = len(AP)
lag = len(AG)



for apidx,ap in enumerate(AP):
    S.tx = RadioNode(typ='tx',name=ap)
    S.tx.loadini(ap+'.ini',rep=pstruc['DIRNETSAVE'])
    for agidx,ag in enumerate(AG):
        print '---------------------'
        print ' Raytracing for :    '
        print ' AP #',AP[apidx-1] ,' / AG #',ag
        print '---------------------'
        print ' Computed :'
        print 'AP:',apidx-1,'/',lap+1
        print 'AG:',agidx,'/',lag
        print '---------------------'
        S.rx = RadioNode(typ='rx',name=ag)
        S.rx.loadini(ag+'.ini',rep=pstruc['DIRNETSAVE'])
        S.run(apidx+1,range(1,S.rx.N+1))

#### STEP 2 : all mobile/mobile

icag = itertools.combinations(AG,2)
for cag in icag:
    S.tx = RadioNode(typ='tx',name=cag[0])
    S.tx.loadini(cag[0]+'.ini',rep=pstruc['DIRNETSAVE'])
    S.rx = RadioNode(typ='tx',name=cag[1])
    S.rx.loadini(cag[1]+'.ini',rep=pstruc['DIRNETSAVE'])
    lidxpts = range(1,S.rx.N+1)
    print '---------------------'
    print ' Raytracing for :    '
    print ' AG #', cag[0] ,' / AG #',cag[1]
    print '---------------------'
    for n in lidxpts:
        print ' trajectory point #',n,'/',S.rx.N+1
        S.run(n,n)



