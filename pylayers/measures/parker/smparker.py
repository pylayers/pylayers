#!/usr/bin/python
#-*- coding:Utf-8 -*-
from serial import Serial
import pdb
import binascii
import matplotlib as plt


class Profile(object):
    accmax = 200
    vmax = 10
    dmax = 1500

    def __init__(self,num,aa,ad,d,v,vs):
        """
        """
        assert(0<aa<accmax)
        assert(0<v<vmax)
        assert(0<d<dmax)

        self.cmd = 'PROFILE'+str(num)+'('+str(aa)+','+str(ad)+',',str(d),+','+str(v)+','+str(vs)

    def show(self):
        """
        """
        pass

class Axes(object):
    def __init__(self,_id,name,scale=12820,typ='t'):
        """
        _id  : axes id
        name : axes name
        scale : nstep/cm if typ='t'

        """
        self.status = {'cpa':0,
                       'loop':0,
                       'wft':0,
                       'runp':0,
                       'me':0,
                       'sta':0
                      }
        self.lprofile=[]

    def add_profile(self,profile):
        """ add a profile in the list
        """
        if len(self.lprofile)<8:
            self.lprofile.append(profile)

    def show(self):
        """
        """
        pass

class StepMotor(object):
    svar  = {'BU':'Buffer Usage',
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
            'UF':'User Program Fault status'
           }


    def __init__(self,port):
        self.ser = Serial(port = port, baudrate=9600, timeout = 1)
        self.axes  = [Axes(1,'x'), Axes(2,'y'), Axes(3,'rot') ] #self.a4  = Axes(4,'z')


    def com2(self,com='1R(ST)'):
        self.ser.write(com+'\r\n')
        st = self.ser.readlines()
        return(st)

    def com(self,axe,name,rg=''):
        if rg!='':
            cst = str(axe)+name+str(rg)+'\r\n'
        else:
            cst = str(axe)+name+'\r\n'

        self.ser.write(cst)
        st = self.ser.readlines()
        return(st)

    def add_profile(axis,aa,ad,d,v,vs):
        ax = self.axis[axis-1]
        npro = len(ax.lprofile)
        prof = Profile(npro+1,aa,ad,d,v,vs)
        ax.add_profile(prof)

    def mvpro(self,axis,id_pro):
        """ move according to specified profile
        Parameters
        ----------

        axis
        id_pro

        Returns
        -------
        str

        Examples
        --------

        sm.move(1,1)  # axis 1 , profile 1
        """

        com = 'USE('+str(id_pro)+')'
        self.com(axis,com)
        self.com(axis,'G')

    def status(self,axis=1):
        st = self.com(axis,'R','(ST)')
        status = st[1]
        status = status.replace('*','').replace('\r\n','').split('_')
        self.axes[axis-1].status['cpa']=status[0][0]
        self.axes[axis-1].status['loop']=status[0][0]
        self.axes[axis-1].status['wft']=status[0][0]
        self.axes[axis-1].status['runp']=status[0][0]
        self.axes[axis-1].status['me']=status[2][0]
        self.axes[axis-1].status['sta']=status[4][3]
        #self.axes[axis-1].status['cpa']=status[0][0]
        #self.axes[axis-1].status['cpa']=status[0][0]
        return status

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
    #st = sm.com(1,'ON')
    #st = sm.com(1,'LIMITS',(0,1,1))
    #st = sm.com(1,'1D-4000')
    #st = sm.com(1,'G')
