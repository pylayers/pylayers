import socket
import time
import struct
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
from types import *
from numpy import array
import pdb
import select
import pylayers.signal.bsignal as bs
from pylayers.signal.channel import *
from time import sleep
"""

Module to drive the network analyzer E5072A

"""

# Adapted from  mclib by Thomas Schmid (http://github.com/tschmid/mclib)

class SCPI:
    PORT = 5025
    _chunk = 128
    _verbose = False
    _timeout = 0.150

    def __init__(self,host,port=PORT,timeout=None,verbose=True):
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

    def _write(self, cmd):
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

    def close(self):
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

    def select(self,param='S11',chan=1):
        """ select parameter on channel chan

        Parameters
        ----------

        param : string
            {'S11','S12','S21','S22'}
        chan : int
            default 1

        Notes
        -----

        working

        """
        com = ":CALC"+str(chan)+":PAR:DEF "+param+"\n"
        self.s.send(com)

    def reset(self):
        """
        Resets the device to known state (with *RST) and clears the error
        log
        """
        # reset and clear device
        #self.s.send(":SYST:FPReset\n")
        #self.s.send("*RST;*CLS", False)
        self.s.send(":SYST:PRES\n")
        #self.s.send(":STOP\n")

    def setntrace(self,chan=1,ntrace=2):
        """ set number of traces
        Parameters
        ----------

        chan : int
        ntrace : 2

        """
        com1 = ":CALC"+str(chan)+":PAR:COUN "+str(ntrace)+" \n"
        self.s.send(com1)

    def setpar(self,chan=1,par='S11'):
        """ set parameter

        Parameters
        ----------

        chan : 1
        par : string
            {'S11','S12','S21','S22'}

        """
        com1 = ":CALC"+str(chan)+":PAR "+par+"\n"
        self.s.send(com1)

    def autoscale(self,win=1,tr=1):
        """ autoscale on window win trace tr
        """
        com ="DISP:WIND"+str(win)+":TRAC"+str(tr)+":Y:SCAL:AUTO"
        self.write(com)

    def setfstar(self,fstartGHz=1.8,sens=1):
        """ frequency start

        Parameters
        ----------

        fstartGHz : float
        fstop : float

        """
        com1 = ":SENS"+str(sens)+":FREQ:START "

        f1 = str(fstartGHz)+"e9\n"

        self.s.send(com1+f1)

        return(fstartGHz)

    def setfstop(self,fstopGHz=2.2,sens=1):
        """ frequency stop

        Parameters
        ----------

        fstartGHz : float
        fstop : float

        """
        com2 = ":SENS"+str(sens)+":FREQ:STOP "

        f2 = str(fstopGHz)+"e9\n"

        self.s.send(com2+f2)

        return(fstopGHz)


    def getnpoints(self,sens=1):
        """ get number of points

        Parameters
        ----------

        sens :

        """
        com1 = ":SENS"+str(sens)+":SWE:POIN?\n"
        self.s.send(com1)
        try:
            npoints = eval(self.s.recv(8).replace("\n",""))
            return npoints
        except socket.timeout:
            print "problem for getteing number of points"

    def setnpoint(self,Npoints=201,sens=1,echo=False):
        """Change the numbers of points

        Parameters
        --------------
        Npoints
        """

        com = "SENS"+str(sens)+":SWE:POIN "+str(Npoints)
        if echo:
            print com
        self.write(com)


    def getfreq(self,sens=1):
        com = ":SENS"+str(sens)+":FREQ:DATA?\n"
        buf = self.read(com)
        f = np.frombuffer(buf,'>f8')
        freq = f[1:]
        return(freq)

        #tab = []
        #while len(tab)<> Npoints:
        #    com2 = ":SENS"+str(sens)+":FREQ:DATA?\n"
        #    self.s.send(com2)
        #    tab = self.s.recv(Npoints*20).split(",")
        #freq = map(lambda x: eval(x),tab)


    def getIdent(self):
        """
        get VNA Identification
        """
        self.s.send("*IDN?\n")
        try:
            data = self.s.recv(1024)
            return data
        except socket.timeout:
            return ""

    def write(self,com):
        self.s.send(com+"\n")

    def ask(self,com):
        com1 = com+"?\n"
        self.s.send(com1)
        try:
            data = self.s.recv(1024)
        except socket.timeout:
            return ""
        return(data)


    def getdata(self,chan=1,Npoints=201):
        """  getdata

        Parameters
        ----------

        chan : int
            channel number
        """


        #self.write("TRIG:SING")
        comm = 'CALC'+str(chan)+':DATA:SDAT?'
        buff = ''
        while len(buff)<>(Npoints*16+8):
            buff = self.read(comm)

        S = np.frombuffer(buff[8:Npoints*16+8],dtype='>f8')
        Y = S.reshape(Npoints,2)
        Y = Y[:,0]+1j*Y[:,1]
        fGHz = self.getfreq()*1e-9
        S21 = bs.FUsignal(x=fGHz,y=Y)
        return S21

if __name__=='__main__':
    vna = SCPI("129.20.33.201",verbose=False)
    ident = vna.getIdent()
    print "Talking to : ",ident
    vna.write("FORM:DATA REAL")
    vna.select(param='S21',chan=1)
    vna.write(":SENS1:SWE:POIN 1201")
    vna.write("DISP:WIND1:TRAC1:Y:SCAL:AUTO")
    vna.s.send(":SENS1:SWE:POIN?\n")
    Npoints = eval(vna.s.recv(56).replace('\n',''))
    print "Npoints : ",Npoints
    # set fmin fmax
    vna.write(":SENS1:FREQ:STAR 1.8e9")
    vna.write(":SENS1:FREQ:STOP 2.2e9")

    #get frequency range
    com = ":SENS1:FREQ:DATA?\n"

    vna.write("TRIG:SING")

    time.sleep(1)
    com1 = ":CALC1:DATA:SDAT?\n"
    #u = np.arange(0,Npoints)*2
    #v = np.arange(0,Npoints)*2+1
    N = 1
    fGHz = np.linspace(1.8,2.2,1201)
    H = FUchannel(fGHz)
    for k in range(N):
        B = vna.read(com1)
        S = np.frombuffer(B[0:Npoints*16],dtype='>f8')
        H.load(S)
        #S21 = S[u]+1j*S[v]
        #try:
        #    res=np.vstack((res,S21.T))
        #except:
        #    res=S21.T
#
#
#    tab = vna.read(com)
#    f = np.frombuffer(tab,'>f8')
#    freq = f[1:]
#    plt.plot(freq)
    vna.close()
    #plt.imshow(abs(res))
    #plt.axis('tight')
    #plt.show()
