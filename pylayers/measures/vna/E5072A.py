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
        # reset and clear device
        self.s.send(":SYST:FPReset\n")

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

    def autoscale(self,win=1):
        pass

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


    def getfreq(self,sens=1,Npoints=201):
        tab = []
        while len(tab)<> Npoints:
            com2 = ":SENS"+str(sens)+":FREQ:DATA?\n"
            self.s.send(com2)
            tab = self.s.recv(Npoints*20).split(",")
        freq = map(lambda x: eval(x),tab)
        return(freq)


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

    def getdata(self,N=1,chan=1,param='S11'):
        """

        """

        TS = []
        self.select(param,chan)
        Npoints = self.getnpoints()
        Npoints=201
        print Npoints
        self.write("FORM:DATA REAL")
        S = np.zeros(2*201,dtype='f2')
        for k in range(N):
            self.write("TRIG:SING")
            self.write('CALC'+str(chan)+':DATA:SDAT?')
            self.s.recv_into(S,4*201*16)
            endS = self.ask("*OPC")

            #S = np.array(map(lambda x: eval(x),S.split(',')))
            #pdb.set_trace()
            S = S.reshape(Npoints,2)
            S = S[:,0]+1j*S[:,1]
            #try:
            #    TS = np.vstack((TS,S))
            #except:
            #    TS = S
        return(S)

#    def measure(self, channel):
#
#        self.s.send(":DIGitize CHANnel" + str(channel) + "\n")
#        self.s.send(":WAVeform:DATA?\n")
#
#        try:
#            data = ""
#        while (not data.endswith('\n')):
#            data += self.s.recv(1024)
#        
#        except socket.timeout:
#            return None
#
#        data = data.split(' ')
#        data = data[1:len(data)]
#	
#
#     	data = np.array(data, dtype='S11')
#
#	    voltage = data.astype(np.float)
#
#        return voltage
#
