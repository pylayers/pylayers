import socket
import doctest
import time
import struct
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
from types import *
from numpy import array
import pdb
import select
from pylayers.util.project import  *
import pylayers.signal.bsignal as bs
import pylayers.antprop.channel as ch
from time import sleep
import seaborn as sns
import os
"""

Module to drive the network analyzer E5072A
Adapted from  mclib by Thomas Schmid (http://github.com/tschmid/mclib)
Enhanced by M.D.BALDE and B.UGUEN

.. currentmodule:: pylayers.measures.vna.E4972A

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

"""

class SCPI(PyLayers):
    PORT = 5025
    _chunk = 128
    #_chunk = 256
    _verbose = False
    _timeout = 0.150

    def __init__(self,port=PORT,timeout=None,verbose=False):
        """
        Parameters
        ----------

        host : ip address
        port : port
        timeout: float
        verbose : boolean

        """
        if "VNA_IP" in os.environ:
            host = os.environ["VNA_IP"]
        else:
            print "VNA IP not defined"
            exit
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
                 print 'SCPI>> connect({:s}:{:d}) failed {:s}',format(host,port,e)
            else:
                raise e

        self.getIdent()
        print self.ident
        assert('E5072A' in self.ident), "E5072A not responding"
        self.points()
        self.freq()
        self.parS()
        self.avrg()
        self.ifband()

    def __repr__(self):
        st = ''
        st = st + '----------------------------'+'\n'
        st = st + '      PARAMETERS            '+'\n'
        st = st + '----------------------------'+'\n'
        st = st + "Talking to         : " + str(self.ident)+'\n'
        st = st + "Channel            : " + str(self.chan) +'\n'
        st = st + "S Parameter        : " + self.param +'\n'
        st = st + "fmin (GHz)         : " + str(self.fGHz[0])+'\n'
        st = st + "fmax (GHz)         : " +str(self.fGHz[-1])+'\n'
        st = st + "Bandwidth (GHz)    : " + str(self.fGHz[-1]-self.fGHz[0])+'\n'
        st = st + "Nbr of freq points : " + str(self.Nf)+'\n'
        st = st + "Avering            : " + self.b +'\n'
        st = st + "Nbr of averages    : " + str(self.navrg)+'\n'
        st = st + "IF Bandwidth (Hz)  : " + str(self.ifbHz)+'\n'
        return(st)

    def _write(self, cmd):
        """ socket write
        """
        if self.s is None: raise IOError('disconnected')

        for i in xrange(0, len(cmd), self._chunk):
            if (i+self._chunk) > len(cmd): idx = slice(i, len(cmd))
            else: idx = slice(i, i+self._chunk)
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

    def ask(self,com):
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
        self.s.close()

    def _read(self):
        if self.s is None: raise IOError('disconnected')
        buf = bytearray()
        data = True
        while data:
            r,w,e = select.select([self.s], [], [self.s], self._timeout)

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

    def parS(self,param='S21',chan=1,tr=1,cmd='get'):
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
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> vna.parS(param='S21',cmd='set')
        >>> vna.close()


        working

        """
        self.param = param
        self.chan  = chan

        #co = ":CALC"+str(chan)+":PAR:DEF"
        co = ":CALC"+str(chan)+":PAR"+str(tr)+":DEF"
        com = co +' '+param
        com1 = com+"\n"

        if cmd == 'get':
            comg = co+'?\n'
            self.s.send(comg)
            #c = self.read(comg)
            #return(c)

        if cmd=='set':
            self.s.send(com1)


    def reset(self):
        """ Resets the device to known state (with *RST) and clears the error
        log

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> vna.close()

        """

        #self.s.send("CALC:PAR:DEL:ALL") #Deletes all measurements on the VNA
        self.s.send(":SYST:PRES\n")

    def trace(self,chan=1,ntrace=1,cmd='get'):
        """ allows get|set the  number of traces.
            traces are a series of measured data points.
            limits of traces : max nbr of win x max nbr of traces per window (24)

        Parameters
        ----------

        chan : int
        ntrace : 2
        cmd    : get|set

        Examples
        --------

        >>> #from pylayers.measures.vna.E5072A import *
        >>> #vna = SCPI()
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> #vna.trace(chan=1,param='S21',ntrace=1,cmd='set')
        >>> #vna.close()


        """
        self.chan  = chan

        #co = ":CALC"+str(chan)+":PAR"+str(tr)+":DEF"
        com = ":CALC"+str(chan)+":PAR:COUN"

        if cmd == 'set':
            coms = com+'  '+str(ntrace)+"\n"
            self.s.send(coms)
        if cmd == 'get':
            comg = com+'?\n'
            self.s.send(comg)

    def autoscale(self,win=1,tr=1):
        """ autoscale on window win trace tr

        Parameters
        ----------

        win : integer
        tr  : integer

        """
        com ="DISP:WIND"+str(win)+":TRAC"+str(tr)+":Y:SCAL:AUTO"
        self.write(com)
    

    def points(self,value=1601,cmd='get',sens=1,echo=False):
        """ get|set  number of points

        Parameters
        ----------

        sens : int
        cmd  : 'get' | 'set'

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> vna.points(201,cmd='set')
        >>> vna.close()


        """
        com = ":SENS"+str(sens)+":SWE:POIN"
        self.Nf = value

        if cmd=='get':
            comg = com+"?\n"
            self.s.send(comg)
            try:
                self.Nf = eval(self.s.recv(8).replace("\n",""))
            except socket.timeout:
                print "problem for getting number of points"

        if cmd=='set':
            coms = com+' '+str(value)
            if echo:
                print coms
            self.write(coms)

    def freq(self,sens=1,fminGHz=1.8,fmaxGHz=2.2,cmd='get'):
        """ get | set frequency ramp

        Parameters
        ----------

        sens : 1
        fminGHz : frequency start (float)
        fmaxGHz : frequency stop  (float)
        cmd     : 'get' | 'set'

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> vna.freq(fminGHz=3.8,fmaxGHz=4.2,cmd='set')
        >>> vna.close()

        """

        if cmd=='get':
            com1 = ":FORM:DATA REAL"
            self.read(com1)
            com2 = ":SENS"+str(sens)+":FREQ:DATA?\n"
            com3 = self.read(com2)
            buf = com3[8:(self.Nf-1)*16+8]
            f = np.frombuffer(buf,'>f8')
            self.fGHz = f/1e9
            fminGHz = self.fGHz[0]
            fmaxGHz = self.fGHz[-1]
            dfGHz = fmaxGHz - fminGHz

        if cmd=='set':
            com1 = ":SENS"+str(sens)+":FREQ:START "
            com2 = ":SENS"+str(sens)+":FREQ:STOP "
            self.fGHz = np.linspace(fminGHz,fmaxGHz,self.Nf)
            f1 = str(fminGHz)+"e9\n"
            f2 = str(fmaxGHz)+"e9\n"

            self.s.send(com1+f1)
            time.sleep(1)
            self.s.send(com2+f2)

    def getIdent(self):
        """ get VNA Identification
        """
        self.s.send("*IDN?\n")
        try:
            #data = self.s.recv(1024)
            self.ident = self.s.recv(1024)
            #return data
        except socket.timeout:
            return ""


    def getdata(self,chan=1,Nmeas=10):
        """ allows getdata from VNA

        Parameters
        ----------
        Nmeas : number of times of measures
        chan  : int
               channel number
        """

        com = 'CALC'+str(chan)+':DATA:SDAT?'
        for k in  range(Nmeas):
            tic = time.time()
            buff = ''
            while len(buff)<>(self.Nf*16+8):
                buff = self.read(com)

            S = np.frombuffer(buff[8:self.Nf*16+8],dtype='>f8')
            Y = S.reshape(self.Nf,2)
            H = Y[:,0]+1j*Y[:,1]
            try:
                tH = np.vstack((tH,H[None,:]))
            except:
                tH = H[None,:]
                #print tH.shape
            #pdb.set_trace()
            S21 = ch.Tchannel(x=self.fGHz,y=tH)
            return S21
            self.close()
            #toc = time.time()
            #t = toc-tic
            #print "Time measurement (ms) :",t



    def avrg(self,sens=1,b='OFF',navrg=16,cmd='getavrg'):
        """ allows get|set the point averaging

        Parameters
        ----------
        b      : boolean (ON/OFF)
        cmd    : getavgr (0 average OFF
                          1 average ON)
                 setavgr
                 getnavgr (preset value = 16)
                 setnavgr
        navgr    : range of average : [1,999]

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> vna.reset()
        >>> vna.freq(fminGHz=2.8,fmaxGHz=3.2,cmd='set')
        >>> vna.avrg()
        >>> vna.avrg(b='ON',cmd='setavrg')
        >>> vna.avrg()
        >>> vna.avrg(navrg=100,cmd='setavrg')


        """
        self.b     = b
        self.navrg = navrg

        co1  = ":SENS"+str(sens)+":AVER"
        co2  = ":SENS"+str(sens)+":AVER:COUN"
        com1 = co1 +' '+b
        com2 = co2 +' '+str(navrg)

        if cmd == 'getavrg':
            com = co1+"?\n"
            self.s.send(com)
            #c = self.read(com)
            #return(c)

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

        ifbHz  : IF Bandwidth (default : 70000Hz)
        cmd     : get|set

        Examples
        --------

        >>> from pylayers.measures.vna.E5072A import *
        >>> vna = SCPI()
        Agilent Technologies,E5072A,MY51100293,A.01.04
        <BLANKLINE>
        >>> vna.ifband(sens=1,ifbHz=70000,cmd='set')
        >>> vna.close()

        """
        self.ifbHz   = ifbHz

        co  = ":SENS"+str(sens)+":BAND"
        com = co +' '+str(ifbHz)

        if cmd == 'get':
            com = co+"?\n"
            self.s.send(com)

        if cmd == 'set':
            com = com+"\n"
            self.s.send(com)


if __name__=='__main__':
    doctest.testmod()
#
#    vna = SCPI(vna_ip,verbose=False)
#    ident = vna.getIdent()
#    #lNpoints = ['201','401','601','801','1601']
#
#    #lNpoints = [1601]
#    print "Talking to : ",ident
#    vna.write("FORM:DATA REAL")
#    #vna.write("SENS:AVER:ON")
#    vna.select(param='S21',chan=1)
#    vna.setf(startGHz=1.8,stopGHz=2.2)
#    lav = [1,999] #average
#    #lsif = ['1000','300000','500000'] #IF band
#    lsif = ['1000'] #IF band
#    lS = []
#    lt = []
#    Npoints = 1601
#    #for Npoints in lNpoints:
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
#        #get frequency range
#        #com = ":SENS1:FREQ:DATA?\n"
#        #tab = vna.read(com)
#        #f = np.frombuffer(tab,'>f8')
#        #freq = f[1:]
#
#    vna.close()
#
#    a0 = np.abs(lS[0].y) # N x Npoints ; IF = 100 KHz
#    #a1 = np.abs(lS[1].y) # N x Npoints ; IF = 300 KHz
#    #a2 = np.abs(lS[2].y) # N x Npoints ; IF = 500 KHz
#    plt.plot(a0[0],label='IF 100KHz')
#    #plt.plot(a1[1],label='IF 300KHz')
#    #plt.plot(a2[2],label='IF 500KHz')
#
#    sns.set_style("darkgrid")
#    plt.xlabel('points')
#    plt.ylabel('Amplitude')
#    plt.title('Evolution of S21 over number of points')
#    plt.legend(loc='best')
#
#    #Variance error
#    #v0=np.var(lS[0].y,axis=0) # N x Npoints
#    #v1=np.var(lS[1].y,axis=0)
#    #plt.semilogy(v0,'b')
#    #plt.semilogy(v1,'r')
#    #sns.tsplot(data=np.abs(S21.y),time=S21.x,err_style="ci_bars")
#    #sns.tsplot(data=np.abs(S21.y),time=S21.x,err_style="ci_band")
#
