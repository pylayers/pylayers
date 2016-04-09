#-*- coding:Utf-8 -*-
import socket
import doctest
import time
import struct
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
from types import *
from numpy import array
import ipdb
import h5py
import select
from pylayers.util.project import  *
import pylayers.signal.bsignal as bs
import pylayers.antprop.channel as ch
from pylayers.util import pyutil as pyu
from pylayers.measures.exploith5 import Mesh5 
import pylayers.measures.switch.ni_usb_6501 as sw
#from  pylayers.measures.parker.smparker import *
from time import sleep
import seaborn as sns
import os
import ConfigParser

"""
Module to drive the network analyzer E5072A
Adapted from  mclib by Thomas Schmid (http://github.com/tschmid/mclib)
Enhanced by M.D.BALDE

.. currentmodule:: pylayers.measures.vna.E5072A

.. autosummary::
    :toctree: generated

SCPI Class
==========

.. autosummary::
    :toctree: generated/

    SCPI.__init__
    SCPI.__repr__
    SCPI._write
    SCPI._read
    SCPI.write
    SCPI.read
    SCPI.ask
    SCPI.close
    SCPI.parS
    SCPI.reset
    SCPI.trace
    SCPI.autoscale
    SCPI.points
    SCPI.freq
    SCPI.getIdent
    SCPI.getdata
    SCPI.avrg
    SCPI.ifband
    SCPI.calibh5

"""

class SCPI(PyLayers):
    PORT = 5025
    _chunk = 128
    _verbose = False
    _timeout = 1

    def __init__(self, port=PORT, timeout=None, verbose=False,Nr=1,Nt=1,emulated = False):
        """
        Parameters
        ----------

        host : ip address
        port : port:1
        timeout: float
        verbose : boolean

        """
        self.emulated = emulated
        if "VNA_IP" in os.environ:
            host = os.environ["VNA_IP"]
        else:
            print "VNA IP not defined"
            exit
        if not self.emulated:
            try:
                self.host = host
                self._verbose = verbose
                self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                if timeout is not None:
                    self.s.settimout(timeout)
                    self._timeout = timeout

                self.s.connect((host, port))
            except socket.error as e:
                if self._verbose:
                    print 'SCPI>> connect({:s}:{:d}) failed {:s}', format(host, port, e)
                else:
                    pass
                    #self.emulated = True


        self.Nf   = 201
        self.getIdent()
        print self.ident
        assert('E5072A' in self.ident), "E5072A not responding"
        self.freq()
        self.points()
        self.parS()
        self.avrg()
        self.ifband()
        self.power_level()

        #self.getdata()

        #initialization of the switch
        if (Nr!=1) and (Nt!=1) and not self.emulated:
            self.switch = sw.get_adapter()
            reattach=False
            if not self.switch:
                raise Exception("No device found")

            #self.switch.device
            self.switch.set_io_mode(0b11111111, 0b11111111, 0b00000000)

    def __repr__(self):
        st = ''
        st = st + '----------------------------'+'\n'
        st = st + '      PARAMETERS            '+'\n'
        st = st + '----------------------------'+'\n'
        st = st + "Talking to         : " + str(self.ident)+'\n'
        st = st + "Channel            : " + str(self.chan) + '\n'
        st = st + "S Parameter        : " + self.param + '\n'
        st = st + "fmin (GHz)         : " + str(self.fGHz[0])+'\n'
        st = st + "fmax (GHz)         : " + str(self.fGHz[-1])+'\n'
        st = st + "Bandwidth (GHz)    : " + str(self.fGHz[-1]-self.fGHz[0])+'\n'
        st = st + "Nbr of freq points : " + str(self.Nf)+'\n'
        st = st + "Avering            : " + self.b + '\n'
        st = st + "Nbr of averages    : " + str(self.navrg)+'\n'
        st = st + "IF Bandwidth (Hz)  : " + str(self.ifbHz)+'\n'
        st = st + "Power level (dBm)  : " + str(self.power)+'\n'
        st = st + "Nbr of measures    : " + str(self.nmeas)+'\n'
        return(st)

    def _write(self, cmd):
        """ socket write
        """
        if self.s is None:
            raise IOError('disconnected')

        for i in xrange(0, len(cmd), self._chunk):
            if (i+self._chunk) > len(cmd):
                idx = slice(i, len(cmd))
            else:
                idx = slice(i, i+self._chunk)
            self.s.sendall(cmd[idx])

        return cmd

    def write(self, cmd):
        try:
            return self._write(cmd + b'\n')
        except IOError as e:
            if self._verbose:
                print 'SCPI>> write({:s}) failed: {:s}'.format(cmd.strip(), e)
            else:
                raise e

#    def write(self,com):
#        self.s.send(com+"\n")

    def ask(self, com):
        com1 = com+"?\n"
        self.s.send(com1)
        try:
            data = self.s.recv(1024)
        except socket.timeout:
            return ""
        return(data)


    def close(self):
        """ close socket
        """
        if not self.emulated:
            self.s.close()

    def _read(self):
        if self.s is None:
            raise IOError('disconnected')
        buf = bytearray()
        data = True
        while data:
            r, w, e = select.select([self.s], [], [self.s], self._timeout)

            if r: # socket readable
                data = self.s.recv(self._chunk)
                if data:
                    buf += data
                else: # Socket readable but there is no data, disconnected.
                    data = False
                    self.close()
            else: # no data in socket
                data = False
        return buf

    def read(self, cmd):
        try:
            cmd = self._write(cmd + '\n')
            ans = self._read()
            if self._verbose:
                print '>> {:s} \n<< {:s} \n'.format(cmd.strip(), ans.strip()),
                if ans == '':
                    cmd = self._write('SYST:ERR?\n')
                    err = self._read()
                    print '>> {:s} \n<< {:s} \n'.format(cmd.strip(), err.strip()),
            return str(ans.strip())
        except IOError as e:
            if self._verbose:
                print 'SCPI>> ask({:s}) failed: {:s}'.format(cmd.strip(), e)
            else:
                raise e

    def parS(self, param='S21', chan=1, tr=1, cmd='get'):
        """ set|get the measurement S parameter of a selected channel

        Parameters
        ----------

        cmd   : 'get'| 'set'
        param : string
            {'S11','S12','S21','S22'}
        chan : int
            default 1

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.parS(param='S21',cmd='set')
        >>> vna.close()


        working

        """
        self.param = param
        self.chan  = chan
        if not self.emulated:
            # co = ":CALC"+str(chan)+":PAR:DEF"
            co = ":CALC"+str(chan)+":PAR"+str(tr)+":DEF"
            com = co + ' '+param
            com1 = com+"\n"

            if cmd == 'get':
                comg = co+'?\n'
                self.s.send(comg)
                # c = self.read(comg)
                # return(c)

            if cmd == 'set':
                self.s.send(com1)


    def reset(self):
        """ Resets the device to known state (with *RST) and clears the error
        log

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.close()

        """

        # self.s.send("CALC:PAR:DEL:ALL") #Deletes all measurements on the VNA
        if not self.emulated:
            self.s.send(":SYST:PRES\n")

    def trace(self, chan=1, ntrace=1, cmd='get'):
        """ allows get|set the  number of traces.
            traces are a series of measured data points.
            limits of traces : max nbr of win x max nbr of traces per window (24)

        Parameters
        ----------

        chan   : int
        ntrace : 2
        cmd    : 'get'|'set'

        Examples
        --------

        >>> #from pylayers.measures.vna.E5072A import *
        >>> #vna = SCPI()
        >>> #vna.trace(chan=1,param='S21',ntrace=1,cmd='set')
        >>> #vna.close()


        """
        self.chan  = chan

        # co = ":CALC"+str(chan)+":PAR"+str(tr)+":DEF"
        com = ":CALC"+str(chan)+":PAR:COUN"

        if cmd == 'set':
            coms = com+'  '+str(ntrace)+"\n"
            self.s.send(coms)
        if cmd == 'get':
            comg = com+'?\n'
            self.s.send(comg)

    def autoscale(self, win=1, tr=1):
        """ autoscale on window win trace tr

        Parameters
        ----------

        win : integer
        tr  : integer

        """
        if not self.emulated:
            com = "DISP:WIND"+str(win)+":TRAC"+str(tr)+":Y:SCAL:AUTO"
            self.write(com)

    def trigger(self, cmd = 'bus'):
        """ configure the trigger that control the acquisition of the measurements

        Parameters
        ----------

        cmd : string
            'bus' : actives the trigger via the BUS
            'go'   : generates a trigger immediately and executes a measurement, regardless of the setting 
                     of the trigger mode.
            'test' : when all of pending operations complete, *OPC returns 1.
        

        """
        if not self.emulated:            
            if cmd == 'bus':
                com = ":TRIG:SOUR BUS"
                self.read(com)

            if cmd == 'go':
                com = ":TRIG:SING"
                self.read(com)

            if cmd == 'test':
                com = "*OPC?"
                self.read(com)


    def points(self, value=1601, cmd='get', sens=1, echo=False):
        """ 'get'|'set'  number of points

        Parameters
        ----------

        sens : int
        cmd  : 'get' | 'set'

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.points(201,cmd='set')
        >>> vna.close()


        """
        com = ":SENS"+str(sens)+":SWE:POIN"
        self.Nf = value
        self.fGHz = np.linspace(self.fGHz[0],self.fGHz[-1],self.Nf)
        if not self.emulated:
            if cmd == 'get':
                comg = com+"?\n"
                self.s.send(comg)
                try:
                    self.Nf = eval(self.s.recv(8).replace("\n", ""))
                except socket.timeout:
                    # print "problem for getting number of points"
                    raise IOError('problem for getting number of points')

            if cmd == 'set':
                coms = com+' '+str(value)
                if echo:
                    print coms
                self.write(coms)

    def power_level(self,chan=1,cmd='get',power=10):
        """ allows 'get' | 'set' the the power level transmitted. 

        Parameters
        ----------

        chan   : number of the channel
        pow    : int power level transmitted
        cmd     : 'get'|'set'

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.power_level(cmd='set',power=10)
        >>> vna.close()

        """
        self.power   = power
        
        if not self.emulated:
            if cmd == 'get':
                self.read("SOUR"+str(chan)+":POW?\n")
                
            if cmd == 'set':
                if power < 18:
                    co  = "SOUR"+str(chan)+":POW"
                    com = co +' '+str(power)
                    # assert(npow > 17)," power level must not exceed 17dBm "
                    com = com+"\n"
                    self.s.send(com)
                else:
                    raise IOError('power level must not exceed 17dBm')

    def freq(self, sens=1, fminGHz=1.8, fmaxGHz=2.2, cmd='get'):
        """ get | set frequency ramp

        Parameters
        ----------

        sens    : 1
        fminGHz : frequency start (float)
        fmaxGHz : frequency stop  (float)
        cmd     : 'get' | 'set'

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.freq(fminGHz=3.8,fmaxGHz=4.2,cmd='set')
        >>> vna.close()

        """
        if not self.emulated:
            if cmd == 'get':
                com1 = ":FORM:DATA REAL"
                self.read(com1)
                com2 = ":SENS"+str(sens)+":FREQ:DATA?\n"
                com3 = self.read(com2)
                buf = com3[8:(self.Nf-1)*16+8]
                f = np.frombuffer(buf, '>f8')
                self.fGHz = f/1e9
                fminGHz = self.fGHz[0]
                fmaxGHz = self.fGHz[-1]
                dfGHz = fmaxGHz - fminGHz

            if cmd == 'set':
                com1 = ":SENS"+str(sens)+":FREQ:START "
                com2 = ":SENS"+str(sens)+":FREQ:STOP "
                self.fGHz = np.linspace(fminGHz, fmaxGHz, self.Nf)
                f1 = str(fminGHz)+"e9\n"
                f2 = str(fmaxGHz)+"e9\n"

                self.s.send(com1+f1)
                time.sleep(1)
                self.s.send(com2+f2)
        else:
            self.fGHz = np.linspace(fminGHz,fmaxGHz,self.Nf)

    def getIdent(self):
        """ get VNA Identification
        """
        if not self.emulated:
            self.s.send("*IDN?\n")
            try:
                # data = self.s.recv(1024)
                self.ident = self.s.recv(1024)
                # return data
            except socket.timeout:
                return ""
        else:
            self.ident = 'E5072A emulated vna'


    def getdata(self, chan=1,Nr=1, Nt=1, Nmeas=10,calibration=False,seed=0):
        """ getdata from VNA

        Parameters
        ----------
        Nmeas   : number of measures
        chan    : int
                  channel number
        Nr      : Number of receiver (used only when calibration ==True)
        Nt      : Number of transmitter (used only when calibration ==True)
        Nmeas   : Number of measurements 
        calibration : boolean (used only for emulation)
            if True the measurement correspond to a dataset for calibration purpose
        seed    : In emulated mode the data are produced from random number. seed is 
            for reproducibility of random experiments.

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> vna = SCPI()
        >>> vna.parS(param='S21',cmd='set')
        >>> S21 = vna.getdata()
        >>> vna.close()

        """

        self.nmeas    = Nmeas
        #ipdb.set_trace()
        #print "Nmeas vna",self.nmeas
        if not self.emulated:
            com = 'CALC'+str(chan)+':DATA:SDAT?'  #This command sets/gets the corrected data array.
            for k in  range(Nmeas):
                buff = ''
                while len(buff) != (self.Nf*16+8):
                    self.trigger(cmd='test')
                    buff = self.read(com)
                    #print"buffer size :",len(buff)
                #ipdb.set_trace()
                S = np.frombuffer(buff[8:self.Nf*16+8], dtype='>f8')
                Y = S.reshape(self.Nf, 2)
                H = Y[:, 0]+1j*Y[:, 1]
                try:
                    tH = np.vstack((tH, H[None,:]))
                except:
                    tH = H[None,:]
            tH = tH[:,None,None,:]
        else:
            if not calibration:
                np.random.seed(seed)
                tau = 300 + 10*np.random.rand(Nr,Nt)
                h = ch.TBchannel()
                h.SalehValenzuela()
                H = h.toFD(fGHz=self.fGHz)
                EH = np.sum(H.y*np.conj(H.y),axis=1)/self.Nf
                stn = np.sqrt(EH)/10.
                N = stn*(np.random.rand(Nmeas,Nr,Nt,self.Nf)+1j*np.random.rand(Nmeas,Nr,Nt,self.Nf))
                tH = H.y[:,None,None,:]+N
                tH = tH*np.exp(-2*1j*np.pi*self.fGHz[None,None,None,:]*tau[None,:,:,None])
            else:
                np.random.seed(seed)
                tau = 300 + 10*np.random.rand(Nr,Nt)
                tH = np.exp(-2*1j*np.pi*self.fGHz[None,None,None,:]*tau[None,:,:,None])
                dtau = 7*np.random.rand()
                tH = tH* np.exp(-2*1j*np.pi*self.fGHz[None,None,None,:]*dtau)

        return tH

    def getchan(self,chan=1,Nmeas=10,fminGHz=1.8,fmaxGHz=2.2):
        """ get a Tchannel from VNA

        Parameters
        ----------
        Nmeas   : number of times of measures
        chan    : int
                  channel number
        fminGHz : frequency start (float)
        fmaxGHz : frequency stop  (float)

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> vna = SCPI()
        >>> vna.parS(param='S21',cmd='set')
        >>> S21 = vna.getdata()
        >>> #plt.plot(np.abs(S21.y)[0])
        >>> vna.close()

        """

        self.nmeas    = Nmeas
        self.fGHz[0]  = fminGHz
        self.fGHz[-1] = fmaxGHz
        f             = np.linspace(fminGHz, fmaxGHz, self.Nf)
        if not self.emulated:
            com = 'CALC'+str(chan)+':DATA:SDAT?'
            # tic = time.time()
            for k in  range(Nmeas):
                buff = ''

                while len(buff) <> (self.Nf*16+8):
                    buff = self.read(com)

                S = np.frombuffer(buff[8:self.Nf*16+8], dtype='>f8')
                Y = S.reshape(self.Nf, 2)
                H = Y[:, 0]+1j*Y[:, 1]
                try:
                    tH = np.vstack((tH, H[None,:]))
                except:
                    tH = H[None,:]
            S21 = ch.Tchannel(x=f, y=tH)
            return S21
            # toc = time.time()
            # t   = toc-tic
            # print "Time measurement (ms) :",t

    def avrg(self,sens=1,b='OFF',navrg=20,cmd='getavrg'):
        """ allows get|set the point averaging

        Parameters
        ----------    def ifband(self,sens=1,ifbHz=70000,cmd='get'):
        b        : boolean (ON/OFF)
        cmd      : getavgr (0 average OFF
                          1 average ON)
                 setavgr
                 getnavgr (preset value = 16)
                 setnavgr
        navgr    : range of average : [1,999]

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.reset()
        >>> vna.freq(fminGHz=2.8,fmaxGHz=3.2,cmd='set')
        >>> vna.avrg()
        >>> vna.avrg(b='ON',cmd='setavrg')
        >>> vna.avrg()
        >>> vna.avrg(navrg=100,cmd='setavrg')
        >>> vna.close()

        """
        self.b     = b
        self.navrg = navrg
        if not self.emulated:
            co1  = ":SENS"+str(sens)+":AVER"
            co2  = ":SENS"+str(sens)+":AVER:COUN"
            com1 = co1 + ' '+b
            com2 = co2 + ' '+str(navrg)
            co3 = ":SENS"+str(sens)+":AVER:CLE"

            if cmd == 'getavrg':
                com = co1+"?\n"
                self.s.send(com)
                # c = self.read(com)
                # return(c)

            if cmd == 'clear':
                com3 = co3
                self.s.send(com3)

            if cmd == 'setavrg':
                com = com1+"\n"
                self.s.send(com)

            if cmd == 'getnavrg':
                com = co2+"?\n"
                self.s.send(com)

            if cmd == 'setnavrg':
                com = com2+"\n"
                self.s.send(com)


    def ifband(self,sens=1,ifbHz=70000,cmd='get'):
        """ allows get|set the IF bandwidth

        Parameters
        ----------

        ifbHz   : IF Bandwidth (default : 70000Hz)
        cmd     : 'get'|'set'

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        >>> vna.ifband(sens=1,ifbHz=70000,cmd='set')
        >>> vna.close()

        """
        self.ifbHz   = ifbHz
        if not self.emulated:
            co  = ":SENS"+str(sens)+":bandwidth"
            com = co + ' '+str(ifbHz)

            if cmd == 'get':
                com = co+"?\n"
                self.s.send(com)

            if cmd == 'set':
                com = com+"\n"
                self.s.send(com)

    def stepped_mode(self,sens=1,mode='STEP'):
        """ activate sweep mode

        Parameters
        ----------

        sens : int
        mode : 'stepped'
        STEPped : source frequency is CONSTANT during measurement of eah displayed point. \
                  More accurate than ANALog. Dwell time can be set in this mode.
        """
        self.mode = mode
        if not self.emulated:
                com  = ":SENS"+str(sens)+":SWE:GEN"+' '+mode
                self.read(com)

    def swept_mode(self,sens=1):
        """ activate sweep mode

        Parameters
        ----------

        sens : int

        ANALog : source frequency is continuously RAMPING during measurement of each \ 
                 displayed point. Faster than STEPped. Sweep time (not dwell time) can be \ 
                 set in this mode.
        """
        if not self.emulated:
                self.read(":SENS"+str(sens)+":SWE:GEN ANAL")



    def calibh5(self,
                 Nr = 4,
                 Nt = 8,
                 _filemesh5 = 'measures',
                 _filecalh5 = 'mcalib',
                 _filecal = 'cal_config.ini',
                 _filevna='vna_config.ini',
                 typ = 'full',
                 gcalm = 1,
                 cables=[],
                 author='',
                 comment='',
                 ):
        """  measure a calibration vector and store in h5 file

        Parameters
        ----------

        Nr : int
        Nt : int
        _filemesh5 : string
            measurement data file prefix
        _filecalh5 : string
            multi antenna calibration data file prefix
        _filecal : string
            variable parameters calibration .ini configuration file name
        _filevna : string
            vna .ini configuration file name
        typ : string
            'full' | 'single'
        gcalm : int
            selected calibration group of _filecalh5 (used only if typ=='single')
        cables : list of strings
        author
        comment

        """

        # set config
        # File from : ~/Pylayers_project/meas

        switch = sw.get_adapter()
        if not switch:
            raise Exception("No device found")
        switch.set_io_mode(0b11111111, 0b11111111, 0b00000000) 
        #set the NI USB mode in order to use


        # store calibration vector in a hdf5 file

        # SISO case
        #     + calibration is saved in the measurement file
        #     + calibration variable parameters are read in _filecal
        #
        # MIMO case
        #     if typ is 'full'
        #       + calibration are saved in the calibration file
        #       + calibration variable parameters are read in _filecal
        #     else
        #       + a single channel calibration is saved in the measurement file
        #       + calibration variable parameters are read in _filecalh5
        #

        Mr = Nr
        Mt = Nt
        # SISO case
        if (Nr==1) and (Nt==1):
            fileh5w = pyu.getlong(_filemesh5, pstruc['DIRMES'])+'.h5'
            self.load_config_vna(_filename=_filevna)
            dcal = self.load_calconfig(_filename=_filecal)
        else:
        # MIMO case
            if typ=='full':
                fileh5w = pyu.getlong(_filecalh5, pstruc['DIRMES'])+'.h5'
                self.load_config_vna(_filename=_filevna)
                dcal = self.load_calconfig(_filename=_filecal)
            if typ=='single':
                Mr = 1
                Mt = 1

                fileh5w = pyu.getlong(_filemesh5, pstruc['DIRMES'])+'.h5'
                # dcal is obtained from _filecalh5 and gcal
                dcal = Mesh5(_filecalh5).get_dcal(gcalm)
                # In a future version it would be better to load from the 
                # parameter which are in the MIMO calibration file instead 
                # of thevna_config
                self.load_config_vna(_filename=_filevna)

        f = h5py.File(fileh5w, "a")
        try:
            ldataset = f.keys()
        except:
            ldataset = []

        lcal =  filter(lambda x : 'cal' in x, ldataset)
        calname = 'cal' + str(len(lcal)+1)
        cal = f.create_group(calname)

        tic = time.time()

        for iT in range(Mt):
            print "connect transmitter :", iT +1
            switch.write_port(0,iT)
            for iR in range(Mr):
                switch.write_port(1,iR)
                print "connect receiver :", iR + 1
                c = ""
                while "g" not in c:
                    c = raw_input("Hit g key ")

                for k in dcal:
                    print "---------------------------------------------------------------------"
                    print "                 Configuration Parameters                            "
                    print   dcal[k]
                    print "---------------------------------------------------------------------"
                    
                    t1 = time.time()
                    for k2 in dcal[k]:
                        if k2=='nf':
                            self.points(dcal[k]['nf'], cmd='set')
                            print "set number of points  :",dcal[k]['nf']
                        if k2=='ifbhz':
                            self.ifband(ifbHz=dcal[k]['ifbhz'], cmd='set')
                            print "set number of ifbHz   :",dcal[k]['ifbhz']
                        if k2=='navrg':
                            self.avrg(b='ON',navrg=dcal[k]['navrg'], cmd='setavrg')
                        if k2=='power':
                            self.power_level(cmd='set',power=dcal[k]['power'])

                    #print "activation of the trigger via the BUS"
                    #get Nmeas calibration vector

                    tic_Dmeas = time.time()
                    print "beginning of the acquisition"
                    time.sleep(2)
                    Dmeas = self.getdata(chan=1,Nmeas=dcal[k]['nmeas'],calibration=True)
                    toc_Dmeas = time.time()
                    print "end of the acquisition (s) :",toc_Dmeas-tic_Dmeas

                    if ((iR==0) and (iT==0)):
                        cal.create_dataset(k, (dcal[k]['nmeas'], Mr, Mt, self.Nf), dtype=np.complex64)
                        cal[k].attrs['fminghz']   = self.fGHz[0]
                        cal[k].attrs['fmaxghz']   = self.fGHz[-1]
                        cal[k].attrs['time']      = time.ctime()
                        cal[k].attrs['author']    = author
                        cal[k].attrs['cables']    = cables
                        cal[k].attrs['comment']   = comment
                        cal[k].attrs['param']     = self.param
                        cal[k].attrs['nt']        = Nt
                        cal[k].attrs['nr']        = Nr
                        cal[k].attrs['nmeas']     = dcal[k]['nmeas']
                        
                        if typ=='single':
                            cal[k].attrs['_filecalh5'] =_filecalh5 
                            cal[k].attrs['gcalm']= gcalm
                            
                    cal[k][:,iR,iT,:] = Dmeas[:,0,0,:]
                    cal[k].attrs['nf']        = self.Nf
                    cal[k].attrs['ifbhz']     = self.ifbHz
                    cal[k].attrs['navrg']     = self.navrg
                    cal[k].attrs['power']     = self.power

                    t2 =  time.time()
                    print "duration set up (s) :",t2-t1

        toc = time.time()

        print "---------------------------------------------------------"
        print "                   END of calibration                    "
        print " measurement time (s)                     :",toc-tic
        print "---------------------------------------------------------"

        f.close()

    def meastime(self):
        """
        """
        lNf = [201,401,801,1601]
        lIF = np.arange(1000,300000,1000)
        #lNf = [201]
        #lIF = [1000,2000,30000]
        #lIF = [1000]
        tval = np.empty((len(lNf),len(lIF)),dtype=np.float64)
        Pb = np.empty((len(lNf),len(lIF)),dtype=np.float64)
        for kf,f in enumerate(lNf):
            print f
            for ki,i in enumerate(lIF):
                print i
                self.points(f,cmd='set')
                self.ifband(ifbHz=i,cmd='set')
                time.sleep(3)
                tic = time.time()
                Dmeas = self.getdata(chan=1,Nmeas=10)
                toc = time.time()
                H  = Dmeas[:,0,0,:]
                Hm = np.mean(H,axis=0)
                Hc = H - Hm
                var = np.mean(np.abs(Hc*np.conj(Hc)),axis=0)
                var2  = np.mean(var)
                print "Time :",toc-tic
                print "Noise Power dB:",10*np.log10(var2)
                tval[kf,ki]=toc-tic
                Pb[kf,ki]= var2
        return(tval,Pb)



    def load_calconfig(self,_filename='cal_config.ini'):
        """ load the whole set of calibration config
        """
        filename = pyu.getlong(_filename, pstruc['DIRMES'])
        cal_conf  = ConfigParser.ConfigParser()
        cal_conf.read(filename)
        #pdb.set_trace()
        sections  = cal_conf.sections()
        di        = {}
        for section in sections:
            di[section] = {}
            options = cal_conf.options(section)
            for option in options:
                # int/float value
                try:
                    di[section][option] = eval(cal_conf.get(section, option))
                # string value
                except:
                    di[section][option] = cal_conf.get(section, option)

        return(di)

    def load_config_vna(self,_filename='vna_config.ini'):
        """ load a vna config file from an .ini file

        Parameters
        ----------

        _filename : string
                   file name extension .ini

        Examples
        --------


        """

        # filename : ~/Pylayers_project/meas
        filename = pyu.getlong(_filename, pstruc['DIRMES'])

        vna_conf  = ConfigParser.ConfigParser()
        vna_conf.read(filename)

        sections  = vna_conf.sections()
        di        = {}
        for section in sections:
            di[section] = {}
            options = vna_conf.options(section)
            for option in options:
                # int/float value
                try:
                    di[section][option] = eval(vna_conf.get(section, option))
                # string value
                except:
                    di[section][option] = vna_conf.get(section, option)


        # be careful no capital word in the sections
        # link between vna and file ini

        # section from  vna_config
        # section : stimulus
    

        self.fminGHz = di['stimulus']['fminghz']
        self.fmaxGHz = di['stimulus']['fmaxghz']
        self.Nf = di['stimulus']['nf']
        self.power = di['stimulus']['power']

        # section : response
        self.param   = di['response']['param']
        self.navrg   = di['response']['navrg']
        self.ifbHz   = di['response']['ifbhz']

        # apply configuration setup
        
        self.freq(fminGHz=self.fminGHz, fmaxGHz=self.fmaxGHz, cmd='set')
        self.points(self.Nf, cmd='set')
        self.power_level(cmd='set',power=self.power)
        self.parS(param=self.param, cmd='set')
        self.ifband(ifbHz=self.ifbHz, cmd='set')
        self.avrg(b='ON',navrg=self.navrg,cmd='setavrg')
        self.autoscale()

        #toc_freq = time.time()
        #print "time measurement for frequency set up (s) :", tic-toc_freq
        #time.sleep(1)
        
        #tic = time.time()
        #self.points(self.Nf, cmd='set')
        #toc_pts = time.time()
        #print "time measurement for points set up (s) :", tic-toc_pts
        #time.sleep(1)
        
        #self.parS(param=self.param, cmd='set')
        #toc_parS = time.time()
        #print "time measurement for parS set up (s) :", tic-toc_parS
        #time.sleep(1)
        
        #self.ifband(ifbHz=self.ifbHz, cmd='set')
        
        
        self.autoscale()


if __name__ == '__main__':
    doctest.testmod()

#    vna = SCPI(vna_ip,verbose=False)
#    ident = vna.getIdent()
# lNpoints = ['201','401','601','801','1601']
#
# lNpoints = [1601]
#    print "Talking to : ",ident
#    vna.write("FORM:DATA REAL")
# vna.write("SENS:AVER:ON")
#    vna.select(param='S21',chan=1)
#    vna.setf(startGHz=1.8,stopGHz=2.2)
# lav = [1,999] #average
# lsif = ['1000','300000','500000'] #IF band
# lsif = ['1000'] #IF band
#    lS = []
#    lt = []
#    Npoints = 1601
# for Npoints in lNpoints:
#    for sif in lsif:
#        vna.point(value=Npoints,cmd='set')
#        vna.write(":SENS1:BAND "+sif)
#        print "Npoints : ",Npoints
#        com1 = ":CALC1:DATA:SDAT?\n"
#        N = 100
#        fGHz = np.linspace(1.8,2.2,Npoints)











#        tic = time.time()
#        for k in range(N):
#            S = vna.getdata(Npoints=Npoints)
#            lt.append(time.time())
#            try:
#                S21.append(S)
#            except:
#                S21=S
#        toc = time.time()
#        print toc-tic
#        lt.append(toc-tic)
#        lS.append(S21)
#        del S21
# get frequency range
# com = ":SENS1:FREQ:DATA?\n"
# tab = vna.read(com)
# f = np.frombuffer(tab,'>f8')
# freq = f[1:]
#
#    vna.close()
#
# a0 = np.abs(lS[0].y) # N x Npoints ; IF = 100 KHz
# a1 = np.abs(lS[1].y) # N x Npoints ; IF = 300 KHz
# a2 = np.abs(lS[2].y) # N x Npoints ; IF = 500 KHz
#    plt.plot(a0[0],label='IF 100KHz')
# plt.plot(a1[1],label='IF 300KHz')
# plt.plot(a2[2],label='IF 500KHz')
#
#    sns.set_style("darkgrid")
#    plt.xlabel('points')
#    plt.ylabel('Amplitude')
#    plt.title('Evolution of S21 over number of points')
#    plt.legend(loc='best')
#
# Variance error
# v0=np.var(lS[0].y,axis=0) # N x Npoints
# v1=np.var(lS[1].y,axis=0)
# plt.semilogy(v0,'b')
# plt.semilogy(v1,'r')
# sns.tsplot(data=np.abs(S21.y),time=S21.x,err_style="ci_bars")
# sns.tsplot(data=np.abs(S21.y),time=S21.x,err_style="ci_band")
#
