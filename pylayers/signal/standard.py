import numpy as np
import ConfigParser
import pylayers.util.pyutil as pyu
from pylayers.util.project import *

class Channel(object):
    """ a radio channel abstraction
    """
    def __init__(self,fcGHz,BMHz,GMHz,attmask,speff,crate):
        """
        Parameters
        ----------

        fcGHz : float
            center frequency
        BMHz : float
            effective bandwith
        GMHz : float
            guard frequency (inter channel gap)
        attmask : float
           required attenuation on band edges
        speff : float
            spectral efficiency (bit/s/Hz)
        crate : float
            coding rate

        """
        self.fcGHz = fcGHz
        self.BMHz  = BMHz
        self.GMHz  = GMHz
        self.attmask  = attmask
        self.speff = speff
        self.crate = crate
        self.f     = np.array([fcGHz-(BMHz+GMHz)/2000.,fcGHz+(BMHz+GMHz)/2000.])

    def __repr__(self):
        """
        """
        st = str(self.fcGHz)+': ['+str(self.f[0])+','+str(self.f[1])+']\n'
        return(st)

    def overlap(self,C):
        """ tests wether 2 channels overlap

        Parameters
        ----------

        C : Channel

        Returns
        -------

        ov : int
            1 if overlaping channel
            0 otherwise

        Notes
        -----

        To be used later for optimization purpose

        """
        if ( ( (self.f[1]>C.f[0]) & (self.f[1] < C.f[1]))
           | ( (self.f[0]<C.f[1]) & (self.f[0] > C.f[0]))):
            return(1)
        else:
            return(0)

    def capacity(self,SNRdB):
        """ calculate capacity

         Parameters
         ----------

         SNRdB : SNR in dB

         Returns
         -------

         C : Channel capacity in Mbit/s (M=1e6)

        """
        C = self.BMHz*np.log(1+10**(SNRdB/10.))/np.log(2)
        return(C)

class Wstandard(object):
    """ Wireless standard class
    """
    def __init__(self,name,Nchannel,modulation):
        self.name = name
        self.Nchannel = Nchannel
        self.modulation = modulation

    def __repr__(self):
        st = self.name+'\n'
        st = st+'-------------------\n'
        for k in range(self.Nchannel):
            st = st + str(k+1) +' '+  self.chan[k+1].__repr__()
        return(st)


    def bandplan(self,fstart,SMHz=5,BMHz=20,GMHz=2,attmask=20,speff=2,crate=0.5):
        """
        """
        self.chan={}
        lfcGHz = fstart+np.arange(self.Nchannel)*SMHz/1000.
        for k,fcGHz in enumerate(lfcGHz):
            self.chan[k+1] = Channel(fcGHz,BMHz,GMHz,attmask,speff,crate)

class AP(object):
    """ Access Point

    Attributes
    ----------

    pos : np.array(,3)
        AP position in 3D
    wstd : Wireless standard
    PtdBm : Transmit Power
    channels : list of used channels
    snsdBm : float
        receiver sensitivity
    nant : int
        number of antennas

    """

    def __init__(self,**kwargs):
        """

        Examples
        --------

        >>> import pylayers.signal.standard as std
        >>> AP1 = AP()
        >>> AP1.load()

        """
        defaults = { 'p' : np.array([0,0,1.2]),
            'name': 'default',
            'wstd': 'ieee80211b',
            'channels':[11],
            'PtdBm':0,
            'sensdBm': -94,
            'nant':1,
        }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.name = kwargs['name']
        self.p = kwargs['p']
        self.wstd = kwargs['wstd']
        self.PtdBm = kwargs['PtdBm']
        self.channels = kwargs['channels']
        self.sensdBm = kwargs['sensdBm']
        self.nant = kwargs['nant']

        self.upfstd(self.wstd)

    def upfstd(self,wstd='ieee80211b'):
        """ update from standard

        Parameters
        ----------

        wstd : string
            wireless standard name

        Notes
        -----

        All defined wireless standard are gathered in one  single file `wstd.ini`

        This function updates `self.s` which contains standard specific parameters

        """

        wstandard = ConfigParser.ConfigParser()
        fp = open(pyu.getlong('wstd.ini',pstruc['DIRSIMUL']))
        wstandard.readfp(fp)
        self.dwstd = dict(wstandard.items(wstd))

        fstart =  eval(self.dwstd['fcghzstart'])
        nchan = eval(self.dwstd['nchan'])
        modulation = self.dwstd['modulation']
        BMHz = eval(self.dwstd['bmhz'])
        GMHz = eval(self.dwstd['gmhz'])
        SMHz = eval(self.dwstd['smhz'])
        attmask = eval(self.dwstd['attmask'])
        self.s  = Wstandard(self.wstd,nchan,modulation)
        self.s.bandplan(fstart,SMHz=SMHz,BMHz=BMHz,GMHz=GMHz,attmask=attmask)
        fp.close()


    def load(self,_fileini='defAP.ini'):
        """ loading an access point from file

        Parameters
        ----------

        _fileini : string
            access point description ini file

        """

        self._fileini  = _fileini
        self.config = ConfigParser.ConfigParser()
        fp = open(pyu.getlong(_fileini,pstruc['DIRSIMUL']))
        self.config.readfp(fp)
        self.ap = dict(self.config.items('ap'))
        self.name = self.ap['name']
        self.p = eval(self.ap['pos'])
        self.wstd = self.ap['wstd']
        self.PtdBm = eval(self.ap['ptdbm'])
        self.channels = eval(self.ap['chan'])
        self.sensdBm = eval(self.ap['snsdbm'])
        self.nant = eval(self.ap['nant'])
        fp.close()
        self.upfstd(self.wstd)



    def __repr__(self):
        st = 'Name : '+str(self.name)+'\n'
        st = st + 'position : '+str(self.p)+'\n'
        st = st+ 'Wireless standard : '+str(self.wstd)+'\n'
        st = st+ 'Trasmit Power (dBm) : '+str(self.PtdBm)+'\n'
        st = st+ 'Selected Channels  : '+str(self.channels)+'   '
        for k in self.channels:
           st = st + self.s.chan[k].__repr__()
        st = st+ 'Sensitivity (dBm) : '+str(self.sensdBm)+'\n'
        st = st+ 'Number of antennas : '+str(self.nant)+'\n'
        return(st)


#Wifi11b = Wstandard('IEEE802.11.b',14,'dsss')
#Wifi11b.bandplan(2.412)
#  IoT (Zigbee alternative
#  MIMO 4x4
#Wifi11ah = Wstandard('IEEE802.11.ah',13,'ofdm')
#"Wifi11ah.bandplan(0.903,SMHz=2,BMHz=1,GMHz=1,attmask=20)
#Bluetooth = Wstandard('Bluetooth',79,'gmsk')
#Bluetooth.bandplan(2.4025,SMHz=1,BMHz=1,GMHz=0,attmask=0)
