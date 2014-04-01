import time
import serial
import thread
import pdb

class LPRS(object):


    def __init__(self,portId=0,baud='U4'):
        self.port = '/dev/ttyUSB'+str(portId)
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

