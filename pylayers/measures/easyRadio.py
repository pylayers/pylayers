#!/usr/bin/python
#-*- coding:Utf-8 -*-
import time
import serial
import thread
import pdb
import os

class LPRS(object):
    """ Class for handling LPRS module


    1) Plug Tx and Rx
    2) check on which serial port they appear dmesg | grep tty
    3) Launch an ipython session for each terminal

    Tx session

    >>> from pylayers.measures.easyradio import *
    >>> Tx=LPRS(0)
    >>> Tx.master()


    Rx session

    >>> from pylayers.measures.easyradio import *
    >>> Tx=LPRS(1)
    >>> Tx.slave()


    """

    def __init__(self,portId=0,baud='U4'):
        self.port = os.path.join('/dev','ttyUSB'+str(portId))
        self.baudrate={'U1':2400,
                  'U2':4800,
                  'U3':9600,
                  'U4':19200,
                  'U5':38400,
                  'U6':31250,
                  'U7':76800,
                  'U8':115200}

        self.power={'P0':-2.5,
                    'P1':-1,
                    'P2':1,
                    'P3':2,
                    'P4':3,
                    'P5':4,
                    'P6':4.5,
                    'P7':5,
                    'P8':6,
                    'P9':7}

        self.bandplan = {'b0':0.8697,
                         'b1':0.902,
                         'b2':0.863 }
        self.bandwidth = {'B0':12.5,
                         'B1':25,
                         'B2':50,
                         'B3':100,
                         'B6':150}

        self.ser = serial.Serial(
        port=self.port,
        baudrate= self.baudrate[baud],
        parity=serial.PARITY_NONE,
        stopbits=serial.STOPBITS_ONE,
        bytesize=serial.EIGHTBITS,
        timeout = 0.5
        )

        self.exitflag = 0

    def __repr__(self):
        st = ''
        c = self.cmd('T3')
        st = st + c +'\n'
        powerId = self.cmd('P?').replace('ER_CMD#','')
        power = self.power[powerId]
        st = st + 'RF Power Output : '+str(power)+' dBm\n'
        ChannelNumber = self.cmd('C?').replace('ER_CMD#','')
        bandplan = self.cmd('b?').replace('ER_CMD#','')
        Bandwidth = self.cmd('B?').replace('ER_CMD#','')
        B  = self.bandwidth[Bandwidth]
        f0 = self.bandplan[bandplan]
        C  = eval(ChannelNumber.replace('C',''))
        st = st + 'Band Plan '+bandplan +'\n'
        st = st + 'ChannelNumber ' + ChannelNumber +'\n'
        st = st + 'Bandwidth ' + str(self.bandwidth[Bandwidth])+ ' kHz \n'
        fc = f0 + C*B/1e6+B/(2e6)
        st = st + 'fc : '+str(fc*1000.)+' MHz\n'

        return(st)


    def cmd(self,cmd='C?'):
        """ Send a command to easyRadio module
        """
        if cmd == 'close':
            self.ser.close()


        prefix = "ER_CMD#"
        command = prefix+cmd

        self.ser.write(command)
        while self.ser.outWaiting()>0:
            pass

        out = ''
        out += self.ser.read(1)
        time.sleep(0.005)
        while self.ser.inWaiting() > 0:
            time.sleep(0.005)
            out += self.ser.read(1)
        # command acknowledgement
        if out != '':
            if command == out:
                time.sleep(0.05)
                self.ser.write("ACK")
                while self.ser.outWaiting()>0:
                    pass
                out = ''
                out += self.ser.read(1)
                time.sleep(0.005)
                if self.ser.inWaiting() <>0:
                    while self.ser.inWaiting() > 0:
                        time.sleep(0.01)
                        out += self.ser.read(1)
        return(out)

    def _listen(self,name,timeout):
        listen_txt = ''
        print "Start listening"
        while True:
            while self.ser.inWaiting()>0:
                time.sleep(0.005)
                listen_txt += self.ser.read(1)
            if listen_txt !='':
                print listen_txt
                listen_txt = ''
            if (self.exitflag == 1):
                print "Listening thread exit"
                self.exitflag = 0
                thread.exit()

    def listen(self):
        thread.start_new_thread( self._listen ,("Listen",2,))
        return(0)

    def beacon(self):
        thread.start_new_thread( self._beacon ,("Beacon",2,))
        return(0)

    def flush(self):
        numbytes = self.ser.inWaiting()
        out = ''
        if numbytes>0:
            out = self.ser.read(numbytes)
            if out != '':
                return(out)


    def master(self,delay=0.005):
        """ master send A and wait for B


        while 1:
            M(send A) -> S
            M(Wait for response B)

        """
        cpt = 0
        tic = time.time()
        while True:
            self.send('A')
            time.sleep(delay)
            resp = ' '
            while (resp[-1]<>'B'):
                resp = ' '
                while self.ser.inWaiting()>0:
                #time.sleep(0.005)
                    resp += self.ser.read(1)
            # wait for ack
                if resp[-1]=='B':
                #print resp[-1]
                    try:
                        rssi = int(self.cmd('T8'),16)
                        print "count ",cpt
                        print "rssi  ",rssi
                    except:
                        pass
                    cpt = cpt+1
                if time.time()-tic>3:
                    tic = time.time()
                    break

    def slave(self,delay=0.005):
        """ slave wait for A and send B

        Parameters
        ----------

        """
        cpt = 0
        tic = time.time()
        while True:
            # In buffer is not empty
            resp = ' '
            while self.ser.inWaiting()>0:
                #time.sleep(0.005)
                resp += self.ser.read(1)
            # wait for ack
            print resp
            if (resp[-1]=='A'):
                try:
                    rssi = int(self.cmd('T8'),16)
                    print cpt,rssi
                except:
                    pass
                self.send('B')
                time.sleep(delay)
                cpt = cpt+1

    def STx(self):
        """ Sounder Tx automate

        Switch Tx antennas
        VNA Measurement

        """
        pass

    def SRx(self):
        """ Sounder Rx automate

        Move scanner
        Switch Rx antennas

        """
        state = 'init'
        while True:
            resp = ' '
            while self.ser.inWaiting()>0:
                #time.sleep(0.005)
                resp += self.ser.read(1)


    def send(self,buf='test'):
        self.ser.write(buf)
        while self.ser.outWaiting()>0:
            pass

    def _beacon(self,name,timeout=2):
        print "Start beaconing"
        while True:
            time.sleep(timeout)
            self.send('Beacon')
            if self.exitflag==1:
                print "Beaconing thread exit"
                self.exitflag = 0
                thread.exit()

