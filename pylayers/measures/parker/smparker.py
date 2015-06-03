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


    dstatus = {}
    dstatus[1]='command processing paused'
    dstatus[2]='looping'
    dstatus[3]='wait for trigger'
    dstatus[4]='running program'
    dstatus[5]='going home'
    dstatus[6]='waiting for delay timeout'
    dstatus[7]='registration in progress'
    dstatus[9]='motor energised'
    dstatus[11]='event trigger active until trigger inputs are reset'
    dstatus[12]='input in LSEL not matching label'
    dstatus[13]='-ve limit seen during last move'
    dstatus[14]='+ve limit seen during last move'
    dstatus[19]='moving'
    dstatus[20]='stationnary'
    dstatus[21]='no registration signal seen in registration window'
    dstatus[22]='cannot stop within the defined regisration distance'

    dusrflt = {}
    dusrflt[1]='value is out of range'
    dusrflt[2]='incorrect command syntax'
    dusrflt[3]='last label already in use'
    dusrflt[4]='label of this name not defined'
    dusrflt[5]='missing Z pulse when homing'
    dusrflt[6]='homing failed no signal detected'
    dusrflt[7]='home signal too narrow'
    dusrflt[8]='drive de-energised'
    dusrflt[9]='cannot related end statement to a label'
    dusrflt[10]='program memory buffer full'
    dusrflt[11]='no more motion profiles avaible'
    dusrflt[12]='no more sequence labels avaible'
    dusrflt[13]='end of travel limit hit'
    dusrflt[14]='still moving'
    dusrflt[15]='deceleration error'
    dusrflt[16]='transmit buffer overflow'
    dusrflt[17]='user program nesting overflow'
    dusrflt[18]='cannot use an undefined profile'
    dusrflt[19]='drive not ready'
    dusrflt[22]='save error'
    dusrflt[23]='command not supported by this product'

    ddrvflt={}
    ddrvflt[1]='Composite Fault'
    ddrvflt[2]='Output stage over curent'
    ddrvflt[3]='Output stage over curent'
    ddrvflt[4]='Output stage over curent'
    ddrvflt[5]='Output stage over curent'
    ddrvflt[6]='Output stage over curent'
    ddrvflt[7]='Output stage over curent'
    ddrvflt[8]='Output stage over curent'

    def __init__(self,_id,name,ser,scale=12820,typ='t'):
        """
        _id  : axes id
        name : axes name
        ser  : serial port 
        scale : nstep/cm if typ='t'
        typ : 't'|'r'

        """
        self.status = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self.usrdflt = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self.drvflt = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self._id = _id
        self.name = name
        # serial port
        self.ser = ser
        # list of profiles
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

    def com2(self,com='R(ST)'):
        self.ser.write(str(self._id)+com+'\r\n')
        st = self.ser.readlines()
        return(st)

    def com(self,name,rg=''):
        if rg!='':
            cst = str(self._id)+name+str(rg)+'\r\n'
        else:
            cst = str(self._id)+name+'\r\n'

        self.ser.write(cst)
        st = self.ser.readlines()
        return(st)

    def home(self):
        if axis!=[]:
            cstr = 'HOME1(-,1,-15,100,2)'
        else:
            cstr = 'HOME1(-,1,-15,100,2)'
        self.com2(cstr)
        self.com2('GH')
        err = self.reg('UF')
        return(err)

    def add_profile(aa,ad,d,v,vs):
        """  add new profile to list
        """
        npro = len(self.lprofile)
        prof = Profile(npro+1,aa,ad,d,v,vs)
        self.lprofile.append(prof)

    def reset(self):
        """ reset axis
        """
        self.com('OFF')
        self.com('ON')

    def mvpro(self,id_pro):
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
        self.com(com)
        self.com('G')

    def read(self):
        pass

    def reg(self,typ='ST'):
        """ read boolean quantities in registers  : ST,UF,DF
        """

        buf = self.com('R','('+typ+')')
        buf = buf[1]
        buf = buf.replace('*','').replace('\r\n','').split('_')

        for k in range(8):
            for l in range(4):
                val = eval(buf[k][l])
                if typ=='ST':
                    self.status[k*4+l] = val
                    if val:
                        print Axes.dstatus[k*4+l+1]
                if typ=='UF':
                    self.usrflt[k*4+l] = val
                    if val:
                        print Axes.dusrdflt[k*4+l+1]
                if typ=='DF':
                    self.drvflt[k*4+l] = val
                    if val:
                        print Axes.drvflt[k*4+l+1]


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

class Scanner(Axes):
    def __init__(self,port):
        self.ser = Serial(port = port, baudrate=9600, timeout = 1)
        self.a  = ['',Axes(1,'x',self.ser),
                       Axes(2,'y',self.ser),
                       Axes(3,'rot',self.ser) ] #self.a4  = Axes(4,'z')





if __name__=="__main__":
    s = Scanner('/dev/ttyUSB0')
    #sm.fromfile('prog1')
    #sm.fromfile('AY')
    #Sc[1].com('ON')
    #st = sm.com(1,'LIMITS',(0,1,1))
    #st = sm.com(1,'1D-4000')
    #st = sm.com(1,'G')
