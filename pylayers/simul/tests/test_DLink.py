from pylayers.simul.link import *
import pylayers.signal.waveform as wvf
import pdb
import warnings
import logging

#logging.basicConfig(filename="testDLink.log",level=logging.DEBUG)
#logging.basicConfig(format='%(asctime)s %(message)s')
#logger = logging.getLogger(__name__)

#warnings.filterwarnings("error")
# set the frequency range
fcGHz=5
WMHz = 3000
Nf = 400
fGHz = np.linspace(fcGHz-WMHz*0.5e-3,fcGHz+WMHz*0.5e-3,Nf)
#fGHz = np.array([fcGHz])
# set the layout
filename = 'defstr.lay'
#L=Layout('CEA.lay',bbuild=1)
logger.info(" Load layout %s",filename)
#filename = 'TC2_METIS.lay'
#filename = 'W2PTIN.lay'
L = Layout(filename)
# set the link
logger.info(" Creates DLink %f,%f,%d",fcGHz,WMHz,Nf)
#
# TODO fix this indoor outdoor boolean issue
#
DL = DLink(L=L, fGHz=fGHz, outdoor=False, applywav=True)
#pdb.set_trace()
#DL.L.indoor = True
#DL.L.indoor=False
DL.L.build()
#DL.ca= 23
#DL.cb=14
#DL.ca= 2
#DL.cb= 5
#DL.b = np.array([761.5,1113,1.2])
DL.Aa=Antenna(typ='Omni')
DL.Ab=Antenna(typ='Omni')
DL.a[2]=1.5
DL.b[2]=1.5
#DL.b=np.array([766,1115,1.8])
tic = time.time()
#DL.eval(verbose=True,force=True,bt=False,cutoff=4,threshold=0.1,ra_vectorized=False)
cutoff = 4
threshold = 0.1
logger.info(" evaluation of DLink cutoff:%d threshold:%f",cutoff,threshold)
DL.eval(verbose=False,force=True,bt=True,cutoff=cutoff,threshold=threshold,ra_vectorized=True,nD=1,diffraction=False)
toc = time.time()
print(toc-tic)
#DL.b=np.array([755,1110,1.5])
#DL.eval(force=['sig','ray','Ct','H'],ra_vectorized=True,diffraction=True)
#dist_a_b = np.sqrt(np.sum((DL.a-DL.b)**2))
#
#if DL.R.los:
#    ak0 = DL.H.ak[0]
#    tk0 = DL.H.tk[0]
#   assert tk0*0.3 == dist_a_b, 'invalid distance'
#    lak0 = 20* np.log10(ak0)
#    Friss= 20*np.log10(2.4)+20*np.log10(dist_a_b) + 32.4
#    assert np.allclose(-lak0,Friss,0.1), 'issue in Friss'

# Point outside
#DL.eval(force=['sig','ray','Ct','H'],ra_vectorized=True,diffraction=True,ra_ceil_H=0)
#DL.eval(force=['sig','ray','Ct','H'],ra_vectorized=True,diffraction=True)
# Angular Spread 
doa = DL.H.doa
dod = DL.H.dod
tau = DL.H.taud
E   = DL.H.energy()[:,0,0]
APSthd=np.sum(dod[:,0]*E)/np.sum(E)
APSphd=np.sum(dod[:,1]*E)/np.sum(E)
APStha=np.sum(doa[:,0]*E)/np.sum(E)
APSpha=np.sum(doa[:,1]*E)/np.sum(E)

