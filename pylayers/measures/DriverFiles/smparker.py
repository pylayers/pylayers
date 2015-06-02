#!/usr/bin/python
#-*- coding:Utf-8 -*-
from serial import Serial
import pdb

class StepMotor(object):
    svar = {'BU':'Buffer Usage',
            'CQ':'Command queuing',
            'DF':'Drive Fault status',
            'EI':'Encoder Input',
            'EO':'Encoder signal output',
            'EP':'Encoder Position',
            'ER':'Feedback encoder resolution',
            'EX':'Coms response Style & echo control',
            'IN':'Inputs',
            'IP':'In position flag',
            'IT':'IN position Time',
            'MC':'Motor current',
            'MR':'Motor resolution',
            'MS':'Motor Standby current',
            'MV':'Moving',
            'PA':'Position Absolute',
            'PE':'Position Error',
            'PI':'Position Incremental',
            'RB':'Ready/Busy flag',
            'RM':'Registration Move',
            'RV':'Revision software',
            'SN':'Serial Number',
            'ST':'Status of indexing',
            'UF':'User Program Fault status',
            
           }
    def __init__(self,port):
        self.ser = Serial(port = port, baudrate=9600, timeout = 0.5)


    def com(self,axe,name,rg=''):
        if rg!='':
            cst = str(axe)+name+str(rg)+'\r\n'
        else:
            cst = str(axe)+name+'\r\n'

        self.ser.write(cst)
        st = self.ser.read(100).replace(cst,'')
        return(st)

    def close(self):
        self.ser.close()

    def fromfile(self,cmdfile,dirfile='./DriverFiles'):
        cmdfile = dirfile+'/'+cmdfile
        fd = open(cmdfile,'r')
        lis = fd.readlines()
        for li in lis:
            li=li.replace("\n","\r\n")
            print li
            self.ser.write(li)
        st = self.ser.read(100)
        #while self.ser.inWaiting()>0:
        #    st += self.ser.read(10)
        #self.ser.close()
        return(st)
    #def translation(self,axis,offset):



if __name__=="__main__":
    sm = StepMotor('/dev/ttyUSB0')
    #sm.fromfile('prog1')
    #sm.fromfile('AY')
    st = sm.com(1,'ON')
    st = sm.com(1,'LIMITS',(0,1,1))
    st = sm.com(1,'1D-8000')
    st = sm.com(1,'G')
