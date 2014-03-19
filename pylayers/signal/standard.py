"""

Band Class
===========

.. autosummary::
    :toctree: generated/

    Band.__init__
    Band.channelize

Channel Class
=============

.. autosummary::
    :toctree: generated/

     Channel.__init__
     Channel.__repr__
     Channel.overlap
     Channel.capacity


Wstandard Class
===============

.. autosummary::
    :toctree: generated/

     Wstandard.__init__
     Wstandard.__repr__
     Wstandard.bandplan

AP Class
========

.. autosummary::
    :toctree: generated/

     AP.__init__
     AP.__repr__
     AP.upfstd
     AP.load

"""
import numpy as np
import json
import ConfigParser
import pylayers.util.pyutil as pyu
from pylayers.util.project import *

class Band(dict):
    """ A Band is a structured portion of the spectrum

    A band is subdivided into channels

    """
    def __init__(self,**kwargs):
        """
        """
        defaults = {
        'zone'  : 'Europe',
        'name'  : 'ISM24',
        'fstart'  : 2.412,
        'fstop'  : 2.472,
        'fstep' : 5,
        'bmhz' : 20 }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self['zone']  = kwargs['zone']
        self['name']  = kwargs['name']
        self['fstart']  = kwargs['fstart']
        self['fstop']  = kwargs['fstop']
        self['fstep'] = kwargs['fstep']
        self['bmhz'] = kwargs['bmhz']
        self.channelize()

    def channelize(self):
        """
        """
        self.chan={}
        self.fcghz = np.arange(self['fstart'],self['fstop'],self['fstep']/1000.)
        for k,fc in enumerate(self.fcghz):
            if (fc>=4) & (fc<5):
                channum = int(np.round((fc-4)*200))
            if (fc>=5) & (fc<6):
                channum = int(np.round((fc-5)*200))
            if fc<4:
                channum = k+1
            self.chan[channum] = Channel(fc,self['bmhz'])

    def select(self,lchan):
        """ select a dictionnary of channel from a list of channels
        """
        dchan = {}
        for k,chan in enumerate(lchan):
           dchan[k] = self[chan]
        return(dchan)

    def load(self,bandname,_fileini='spectrum.ini'):
        """ load spectrum
        """
        self._fileini = _fileini
        self.config = ConfigParser.ConfigParser()
        fp = open(pyu.getlong(_fileini,pstruc['DIRSIMUL']))
        self.config.readfp(fp)
        band = dict(self.config.items(bandname))

        self['name'] = bandname
        self['zone'] = band['zone']
        self['fstart'] = eval(band['fstart'])
        self['fstop'] = eval(band['fstop'])
        self['fstep'] = eval(band['fstep'])
        self['bmhz'] = eval(band['bmhz'])
        self.channelize()

class Channel(dict):
    """ a radio channel abstraction
    """
    def __init__(self,fcghz,bmhz,gmhz=0):
        """
        Parameters
        ----------

        fcghz : float
            center frequency
        bmhz : float
            effective bandwith
        gmhz : float
            guard frequency (inter channel gap)

        """
        self['fcGHz'] = fcghz
        self['BMHz']  = bmhz
        self['GMHz']  = gmhz
        self.fghz = np.array([fcghz-(bmhz+gmhz)/2000.,fcghz+(bmhz+gmhz)/2000.])

    def __repr__(self):
        """ representation
        """
        st = str(self['fcGHz'])+': ['+str(self.fghz[0])+','+str(self.fghz[1])+']\n'
        return(st)

    def __add__(self,chan):
        """ add two adjascent channels
        """
        if (self.fghz[1]==chan.fghz[0]):
            self['fcGHz'] = self.fghz[1]
            self['BMHz'] = self['BMHZ']+chan['BMHz']
            self.fghz[1]=chan.fghz[1]
        elif (self.fghz[0]==chan.fghz[1]):
            self['fcGHz'] = self.fghz[0]
            self['BMHz'] = self['BMHZ']+chan['BMHz']
            self.fghz[0]=chan.fghz[0]
        else:
            pass
        return(self)

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
        """ calculates channel capacity

         Parameters
         ----------

         SNRdB : SNR in dB

         Returns
         -------

         C : Channel capacity in Mbit/s (M=1e6)

        """
        C = self.bmhz*np.log(1+10**(SNRdB/10.))/np.log(2)
        return(C)

class Wstandard(dict):
    """ Wireless standard class
    """
    def __init__(self,stdname,_filejson='wstd.json'):
        """
        Parameters
        ----------

        name : string

        Examples
        ---------

        >>> from pylayers.signal.standard import *
        >>> Wifiag =Wstandard('ieee80211ah')

        """
        self.name = stdname
        self.load(stdname)

    def __repr__(self):
        st = self.name+'\n'
        st = st+'-------------------\n'
        for k in range(self.Nchannel):
            st = st + str(k+1) +' '+  self.chan[k+1].__repr__()
        return(st)


    def load(self,stdname,_fileini='wstd.json'):
        """ load a standard from file

        Parameters
        ----------

        stdname : string
            standard name
        _fileini : string
            file containing the description of available standards



        """

        fp = open(pyu.getlong('wstd.json',pstruc['DIRSIMUL']))
        stds = json.load(fp)
        fp.close()
        std = stds[stdname]
        for k in std:
            if k<> "channels":
                try:
                    self[k] = eval(std[k])
                except:
                    self[k] = std[k]
            else:
                chan = std[k]
                for k in chan:
                    bandname = k
                    fstart=chan[k]['fstart']
                    fstop =chan[k]['fstop']
                    smhz =chan[k]['smhz']
                    bmhz = chan[k]['bmhz']
                    gmhz = chan[k]['gmhz']
                    self.bandplan(fstart=fstart,fstop=fstop,smhz=smhz,bmhz=bmhz,gmhz=gmhz)



    def bandplan(self,fstart,fstop,smhz=5,bmhz=20,gmhz=2):
        """ construct the different channels of the standard

        Parameters
        ----------

        fstart : start frequency GHz
        fstop : stop frequency GHz
        smhz : step between adjacscent channels
        bmhz : useful channel bandwidth
        gmhz : gap between channels

        """
        if not(hasattr(self,'chan')):
            self.chan={}
        Nchannel = np.round((fstop-fstart)/(smhz/1000.)).astype(int)+1
        fcghz = np.linspace(fstart,fstop,Nchannel,endpoint=True)
        for k,fc in enumerate(fcghz):
            if (fc>=4) & (fc<5):
                channum = int((fc-4)*200)
            if (fc>=5) & (fc<6):
                channum = int((fc-5)*200)
            if fc<4:
                channum = k+1
            self.chan[channum] = Channel(fc,bmhz,gmhz)
        try:
            self.fcghz=np.hstack((self.fcghz,fcghz))
        except:
            self.fcghz=fcghz

        self.fcghz=np.sort(self.fcghz)


class AP(dict):
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

        self['name'] = kwargs['name']
        self['p'] = kwargs['p']
        self['wstd'] = kwargs['wstd']
        self['PtdBm'] = kwargs['PtdBm']
        self['channels'] = kwargs['channels']
        self['sensdBm'] = kwargs['sensdBm']
        self['nant'] = kwargs['nant']

        self.upfstd(self['wstd'])

    def __repr__(self):
        """ specific representation

        It respects keys of the dictionnary
        """
        st = 'name : '+str(self['name'])+'\n'
        st = st + 'p : '+str(self['p'])+'\n'
        st = st+ 'wstd : '+str(self['wstd'])+'\n'
        st = st+ 'PtdBm : '+str(self['PtdBm'])+'\n'
        st = st+ 'channels  : '+str(self['channels'])+'   '
        for k in self['channels']:
           st = st + self.s.chan[k].__repr__()
        st = st+ 'sensdBm : '+str(self['sensdBm'])+'\n'
        st = st+ 'nant : '+str(self['nant'])+'\n'
        return(st)


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
        fp.close()

        self.dwstd = dict(wstandard.items(wstd))
        fstart =  eval(self.dwstd['fcghzstart'])
        nchan = eval(self.dwstd['nchan'])
        modulation = self.dwstd['modulation']
        BMHz = eval(self.dwstd['bmhz'])
        GMHz = eval(self.dwstd['gmhz'])
        SMHz = eval(self.dwstd['smhz'])
        attmask = eval(self.dwstd['attmask'])
        self.s = Wstandard(wstd,nchan,modulation)
        self.s.bandplan(fstart,SMHz=SMHz,BMHz=BMHz,GMHz=GMHz,attmask=attmask)


        
    def load(self,_fileini='defAP.ini'):
        """ loading an access point from file

        Parameters
        ----------

        _fileini : string
            access point description ini file

        """

        self._fileini = _fileini
        self.config = ConfigParser.ConfigParser()
        fp = open(pyu.getlong(_fileini,pstruc['DIRSIMUL']))
        self.config.readfp(fp)

        ap = dict(self.config.items('ap'))
        self['name'] = ap['name']
        self['p'] = eval(ap['pos'])
        self['wstd'] = ap['wstd']
        self['PtdBm'] = eval(ap['ptdbm'])
        self['channels'] = eval(ap['chan'])
        self['sensdBm'] = eval(ap['snsdbm'])
        self['nant'] = eval(ap['nant'])

        fp.close()
        self.upfstd(self['wstd'])




#Wifi11b = Wstandard('IEEE802.11.b',14,'dsss')
#Wifi11b.bandplan(2.412)
#  IoT (Zigbee alternative
#  MIMO 4x4
#Wifi11ah = Wstandard('IEEE802.11.ah',13,'ofdm')
#"Wifi11ah.bandplan(0.903,SMHz=2,BMHz=1,GMHz=1,attmask=20)
#Bluetooth = Wstandard('Bluetooth',79,'gmsk')
#Bluetooth.bandplan(2.4025,SMHz=1,BMHz=1,GMHz=0,attmask=0)
