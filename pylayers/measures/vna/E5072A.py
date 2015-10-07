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
import pylayers.antprop.channel as ch
from time import sleep
import seaborn as sns
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
        """
        Parameters
        ----------

        host : ip address
        port : port
        timeout: float
        verbose : boolean

        """
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

    def setf(self,startGHz=1.8,stopGHz=2.2,sens=1):
        """ frequency start

        Parameters
        ----------

        startGHz : float
        stopGHz : float

        """
        com1 = ":SENS"+str(sens)+":FREQ:START "
        com2 = ":SENS"+str(sens)+":FREQ:STOP "

        f1 = str(startGHz)+"e9\n"
        f2 = str(stopGHz)+"e9\n"

        self.s.send(com1+f1)
        time.sleep(1)
        self.s.send(com2+f2)



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
            print "problem for getting number of points"

    def setnpoint(self,Npoints=1601,sens=1,echo=False):
        """Change the number of points

        Parameters
        --------------
        Npoints
        """

        com = "SENS"+str(sens)+":SWE:POIN "+str(Npoints)
        if echo:
            print com
        self.write(com)


    def getfreq(self,sens=1):
        """
        Returns
        -------

        fGHz : np.array
        """
        com = ":SENS"+str(sens)+":FREQ:DATA?\n"
        buf = self.read(com)
        f = np.frombuffer(buf,'>f8')
        freq = f[1:]
        fGHz = freq/1e9
        return(fGHz)

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


    def getdata(self,chan=1,Npoints=1601):
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
        fGHz = self.getfreq()
        S21 = ch.Tchannel(x=fGHz,y=Y)
        return S21

if __name__=='__main__':
    vna = SCPI("129.20.33.201",verbose=False)
    ident = vna.getIdent()
    #lNpoints = ['201','401','601','801','1601']

    #lNpoints = [1601]
    print "Talking to : ",ident
    vna.write("FORM:DATA REAL")
    #vna.write("SENS:AVER:ON")
    vna.select(param='S21',chan=1)
    vna.setf(startGHz=1.8,stopGHz=2.2)
    lav = [1,999] #average
    #lsif = ['1000','300000','500000'] #IF band
    lsif = ['1000'] #IF band
    lS = []
    lt = []
    Npoints = 1601
    #for Npoints in lNpoints:
    for sif in lsif:
        vna.setnpoint(Npoints=Npoints)
        vna.write(":SENS1:BAND "+sif)
        print "Npoints : ",Npoints
        com1 = ":CALC1:DATA:SDAT?\n"
        N = 100
        fGHz = np.linspace(1.8,2.2,Npoints)
        tic = time.time()
        for k in range(N):
            S = vna.getdata(Npoints=Npoints)
            lt.append(time.time())
            try:
                S21.append(S)
            except:
                S21=S
        toc = time.time()
        print toc-tic
        lt.append(toc-tic)
        lS.append(S21)
        del S21
        #get frequency range
        #com = ":SENS1:FREQ:DATA?\n"
        #tab = vna.read(com)
        #f = np.frombuffer(tab,'>f8')
        #freq = f[1:]

    vna.close()

    a0 = np.abs(lS[0].y) # N x Npoints ; IF = 100 KHz
    #a1 = np.abs(lS[1].y) # N x Npoints ; IF = 300 KHz
    #a2 = np.abs(lS[2].y) # N x Npoints ; IF = 500 KHz
    plt.plot(a0[0],label='IF 100KHz')
    #plt.plot(a1[1],label='IF 300KHz')
    #plt.plot(a2[2],label='IF 500KHz')

    sns.set_style("darkgrid")
    plt.xlabel('points')
    plt.ylabel('Amplitude')
    plt.title('Evolution of S21 over number of points')
    plt.legend(loc='best')

    #Variance error
    #v0=np.var(lS[0].y,axis=0) # N x Npoints
    #v1=np.var(lS[1].y,axis=0)
    #plt.semilogy(v0,'b')
    #plt.semilogy(v1,'r')
    #sns.tsplot(data=np.abs(S21.y),time=S21.x,err_style="ci_bars")
    #sns.tsplot(data=np.abs(S21.y),time=S21.x,err_style="ci_band")

