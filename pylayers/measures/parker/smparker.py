#!/usr/bin/python
#-*- coding:Utf-8 -*-
from serial import Serial
import pdb
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pylayers.util.project import *

def gettty():
    """get tty and handles port conflicts

    Examples
    --------
    >>> import  os
    >>> from pylayers.measures.parker import smparker
    >>> port = getty()

    """
    import  os
    line = os.popen('dmesg | grep tty | tail -1').read() .replace('\n','')
    tty = line.split('ttyUSB')
    num = tty[1]
    port = '/dev/ttyUSB'+num
    return port

class Profile(PyLayers):
    Accmax = 200
    Vmax = 10
    Dmax = 1500

    def __init__(self,**kwargs):
        """
        Parameters
        ----------

        num     : profile number
        aa      : aceleration (rps2)
        ad      : deceleration (rps2)
        dstep   : distance (step)
        vmax    : vmax (rps)
        vs      : vstart (rps)
        spr     : steps per round


        """

        defaults = {'num': 1,
                    'aa': 200,
                    'ad': 200,
                    'dstep': 12800,
                    'vmax': 15,
                    'vs': 0,
                    'spr': 4000,
                    'N':100}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.num = kwargs['num']
        self.aa = kwargs['aa']
        self.ad = kwargs['ad']
        self.dstep = kwargs['dstep']
        self.vmax = kwargs['vmax']
        self.vs = kwargs['vs']
        self.spr = kwargs['spr']
        self.N = kwargs['N']

        #
        # spr   : steps per revolution
        # dstep : number of steps
        # drev  : number of revolutions
        #

        self.drev = self.dstep/(1.0*self.spr)
        #self.T = (self.drev+self.vmax**2/self.aa)/(1.0*self.vmax)
        self.T = (self.drev+0.5*self.vmax**2/self.aa+0.5*self.vmax**2/self.ad)/(1.0*self.vmax)

        #assert(0<accmax)
        #assert(0<v<vmax)
        #assert(0<d<dmax)

        self.cmd = 'PROFILE'+str(self.num)+\
                '('+str(self.aa)+','+str(self.ad)+\
                ','+str(self.dstep)+','+str(self.vmax)+\
                ','+str(self.vs)+')'

        self.t = np.linspace(0,self.T,self.N)
        self.v = self.vmax*np.ones(self.N)
        self.v1 = self.t*self.aa
        self.v2 = self.t*self.ad

        #
        # t1 : end of acceleration phase
        # t2 : begining of deceleration phase
        #
        t1 = self.vmax/(1.0*self.aa)
        t2 = self.T-self.vmax/(1.0*self.ad)

        u1 = np.where(self.t<t1)[0]
        u2 = np.where(self.t>=t2)[0]

        self.v[u1] = self.t[u1]*self.aa
        self.v[u2] = -self.t[u2]*self.ad+(self.vmax+self.ad*t2)
        #
        # dr : distance in number of revolutions
        # ds : distance in number of steps
        #
        self.dr = np.cumsum(self.v)*(self.t[1]-self.t[0])
        self.ds = self.dr*self.spr

    def __repr__(self):
        st = ''
        st = st + 'num : '+ str(self.num) + '\n'
        st = st + 'acceleration  : '+ str(self.aa) + '\n'
        st = st + 'deceleration  : '+ str(self.ad) + '\n'
        st = st + 'dstep : '+ str(self.dstep) + '/'+ str(self.ds[-1]) + '\n'
        st = st + 'maximum velocity : '+ str(self.vmax) + '\n'
        st = st + 'starting velocity : '+ str(self.vs) + '\n'
        st = st + 'spr : '+ str(self.spr) + '\n'
        st = st + 'N  : '+ str(self.N) + '\n'
        return(st)

    def duration(self):
        """
        """
        self.duration

    def show(self):
        """Enables view profile

        Examples
        --------

        P1=Profile(1,200,200,15000,15,0)
        P1.show()

        """
        plt.subplot(211)
        plt.plot(self.t,self.v)
        plt.xlabel('time (s)')
        plt.ylabel('Velocity (rev/s)')
        plt.title('Evolution of velocity w.r.t time')
        plt.subplot(212)
        plt.plot(self.t,self.ds)
        plt.xlabel('time (s)')
        plt.ylabel('distance (step)')
        plt.title('Evolution of distance w.r.t time')
        plt.legend()
        plt.show()

class Axes(PyLayers):
    """This class allows the control of axes


    svar    : system variables
    dstatus : status bits
    dusrflt : user faults
    ddrvflt : drive faults


    """
    dvar  = {'BU':'Buffer usage',
            'CQ':'Command queuing',
            'DF':'Drive fault status',
            'EI':'Encoder input',
            'EO':'Encoder signal output',
            'EP':'Encoder position',
            'ER':'Feedback encoder resolution',
            'EX':'Coms response style & echo control',
            'IN':'Inputs',
            'IP':'In position flag',
            'IT':'In position time',
            'MC':'Motor current',
            'MR':'Motor resolution',
            'MS':'Motor standby current',
            'MV':'Moving',
            'PA':'Position absolute',
            'PE':'Position error',
            'PI':'Position incremental',
            'RB':'Ready/Busy flag',
            'RM':'Registration Move',
            'RV':'Revision software',
            'SN':'Serial number',
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
    dstatus[13]='-ve : Fin de course atteint'
    dstatus[14]='+ve Fin de course atteint'
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
    ddrvflt[3]='Supply rail failure'
    ddrvflt[4]='Ambient over temperature'
    ddrvflt[5]='Drive over temperature'
    ddrvflt[6]='Configuration error'
    ddrvflt[7]='Motor high voltage rail failure'
    ddrvflt[8]='Output fault'

    def __init__(self,
                 _id=1,
                 name='x',
                 ser=Serial(port=gettty(),baudrate=9600,timeout=0.05),
                 scale=12800,
                 typ ='t'):
        """

        Parameters
        ----------
        _id  : axes id
        name : axes name
        ser  : serial port
        scale : nstep if typ='t'
        typ : 't' | 'r'


        Examples
        --------
        >>> X = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> Y = Axes(2,'y',typ='t',scale=22800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> R = Axes(3,'r',typ='r',scale=2111.111111111111,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        """
        self.status = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self.usrflt = [0,0,0,0,
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
        #
        # serial port
        #
        self.ser = ser
        #
        # list of profiles
        #
        self.scale = scale
        self.typ = typ
        self.lprofile=[]
        self.add_profile()

    def info(self):
        """gives informations about system variables of PARKER

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.info()
        """
        for k in self.dvar:
            st = self.com('R('+k+')')
            print self.dvar[k],st[1].replace('\n','')

    def __repr__(self):

        st = ''
        st = st + '---------------------\n'
        st = st + ' Parameters ' + '\n'
        st = st + '---------------------\n'
        st = st + 'axes_id : ' + str(self._id) + '\n'
        st = st + 'scale : ' + str(self.scale) + '\n'
        st = st + 'typ : ' + str(self.typ) + '\n'

        for k in enumerate(str(self._id)):
            st = st + '---------------------\n'
            st = st + ' Status '+' '+'Axe '+str(k[0]+1)+ '\n'
            st = st + '---------------------\n'
            st1 = self.reg()
            st =  st +st1
        for k,p in enumerate(self.lprofile):
            st = st + '--------------------\n'
            st = st + ' Profile '+str(k+1)+ '\n'
            st = st + '--------------------\n'
            st  = st +  p.__repr__()

        return(st)

    def show(self):
        """
        """
        pass

    #def com2(self,com='R(ST)'):
        #self.ser.write(str(self._id)+com+'\r\n')
        #st = self.ser.readlines()
        #return(st)

    def getvar(self,lvar=[]):
        """allows get state of variables

        Parameters
        ----------

        lvar : list of variables

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.getvar('PA') #Get Position absolute

        """
        if lvar == []:
            lvar = Axes.dvar.keys()
        else:
            var = lvar
            st = self.com('R('+var+')')
            print Axes.dvar[var],st[1]


    #def com(self,command,arg='',verbose=False):
        #""" Send command to serial port

        #Parameters
        #----------

        #command : str command prefix
        #arg :  command argument
        #verbose :

        #"""
        #if arg!='':
            #cst = str(self._id)+command+str(arg)+'\r\n'
        #else:
            #cst = str(self._id)+command+'\r\n'
            #print cst
        #self.ser.write(cst)
        #st = self.ser.readlines()
        #return(st)

    def com(self,command='R(SN)',verbose=False):
        """ send command to serial port

        Parameters
        ----------

        prefix : str command prefix (command='R(SN)')
        arg :  command argument
        verbose :

        Examples
        --------
        >>> # Some suitable commands
        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.com('ON') #Energized motor
        >>> A.com('LIMITS(1,0,0)') #Disable limit +,limits normally open, mode 0 stop motion and abort prog
        >>> A.com('V15') # Change velocity to 15
        >>> A.com('AA15') # Change acceleration to 15
        >>> A.com('AD15') # Change deceleration to 15


        """

        #cst = str(self._id) + command + '\r\n'
        cst = str(self._id) + command + '\r\n'
        self.ser.write(cst)
        st = self.ser.readlines()
        if verbose:
            print cst
        return(st)


    def limits(self,cmd='get',**kwargs):
        """give state and set up limits

        Parameters
        ----------

        cmd  :  'get', 'set'

        """

        defaults = {'mask':0,
                    'typ':1,
                    'mode':1,
                    'LD':200
        }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        mask = kwargs['mask']
        typ  = kwargs['typ']
        mode = kwargs['mode']
        LD   = kwargs['LD']

        if cmd=='get':
            st = self.com('LIMITS')
            ans = st[1].split(' ')

            if '0' in ans[0]:
                print "Enable limits (default setting), "
            if '1' in ans[0]:
                print "Disable limit +, "
            if '2' in ans[0]:
                print "Disable limit -, "
            if '3' in ans[0]:
                print "Disable limit + & -, "

            if '0' in ans[1]:
                print "Limits normally closed (default setting), "
            else:
                print "Limits normally open, "

            if '0' in ans[2]:
                print "Stop motion when a limit is hit and abort the program (default setting), "
            else:
                print "Stop motion when a limit is hit but continue the program, "

            print 'deceleration : ',eval(ans[3].split('D')[1]), "rps²"


        if cmd=='set':
            cstr = 'LIMITS'+'('+str(mask)+','+str(typ)+','+str(mode)+','+str(LD)+')'
            self.com(cstr,verbose=True)

    def home(self,cmd='get',**kwargs):
        """ enables back material home for each axe

        Parameters
        ----------

        cmd  : 'get','set','go'

        Examples
        --------

         >>> s.a[1].home() or A.home() #get informations about the status of HOME
         >>> s.a[1].home('set',vel=15,acc=200) or A.home('set',vel=15,acc=200) #set vel & acc
         >>> s.a[1].home('go') or A.home('go') #back Home  (material)

        """

        defaults = {'mode':0,
                    'vel':15,
                    'acc':20,
                    'edg':'+',
                    'typ':0,
                    'armed':1
                }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.mode = kwargs['mode']
        self.vel = kwargs['vel']
        self.acc = kwargs['acc']
        self.edg = kwargs['edg']
        self.typ = kwargs['typ']
        self.armed = kwargs['armed']

        if cmd=='get':
            st = self.com('HOME')
            ans = st[1].split(' ')

            if '1' in ans[0]:
                print "armed, "
            else:
                print "not armed, "
            if '-' in ans[1]:
                print "reference edge is negative, "
            else:
                print "reference edge is positive, "
            if '1' in ans[2]:
                print "home switch normally closed 1, "
            else:
                print "home switch normally open 0 (default),  "

            if '+' in ans[3]:
                print 'velocity : +',eval(ans[3].split('V+')[1]), "rps"
            else:
                print 'velocity : -',eval(ans[3].split('V-')[1]), "rps"

            print 'acceleration : ',eval(ans[4].split('A')[1]), "rps²"

            if '0' in ans[5]:
                print 'Mode 0: Motor in the active window of the switch(default)'
            if '1' in ans[5]:
                print 'Mode 1: Motor in the position to the edge + or -'
            if '2' in ans[5]:
                print 'Mode 2: Improve homing repeatability'

        if cmd=='set':

            if self.vel>0:
                vel = '+'+str(self.vel)
            else:
                vel = '-'+str(self.vel)


            cstr = 'HOME'+str(self.armed)+\
                          '('+self.edg+','+\
                          str(self.typ)+','+\
                          vel+','+\
                          str(self.acc)+','+\
                          str(self.mode)+')'
            self.com(cstr)

        if cmd=='go':
            cstr = 'HOME1'
            self.com(cstr)
            cstr = 'ARM1'
            self.com(cstr)
            cstr = 'GH'
            self.com(cstr)

        #if lvar == []:
            #lvar = Axes.svar.keys()
        #else:
            #var = lvar['PA']
        #if cmd=='go':
            #if lvar <0:
                #cstr = 'HOME1'
                #self.com(cstr)
                #cstr = 'ARM1'
                #self.com(cstr)
                #cstr = 'GH'
                #self.com(cstr)
            #else:
                #cstr = 'HOME1'
                #self.com(cstr)
                #cstr = 'ARM1'
                #self.com(cstr)
                #cstr = 'H-'
                #self.com(cstr)
                #cstr = 'G'
                #self.com(cstr)


    def add_profile(self,**kwargs):
        """ add a new profile to lprofile

        Parameters
        ----------

        num     : profile number
        aa      : aceleration (rps2)
        ad      : deceleration (rps2)
        dstep   : distance (step)
        v       : vmax (rps)
        vs      : vstart (rps)
        spr     : steps per round

        Examples
        --------


        """
        defaults = {'num': 0,
                    'aa': 200,
                    'ad': 100,
                    'dstep': 12800,
                    'vmax': 15,
                    'vs': 0,
                    'spr': 4000,
                    'N':100}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.num = kwargs['num']
        self.aa = kwargs['aa']
        self.ad = kwargs['ad']
        self.dstep = kwargs['dstep']
        self.vmax = kwargs['vmax']
        self.vs = kwargs['vs']
        self.spr = kwargs['spr']
        self.N = kwargs['N']

        num = len(self.lprofile)
        kwargs['num']=num+1
        if len(self.lprofile)<8:
            prof = Profile(**kwargs)
            # update profile
            self.lprofile.append(prof)
            # send command
            # self.com(prof.cmd)

    def set_profile(self,num):
        """
        """

        assert(num<=len(self.lprofile)),"profile number not defined"
        #self.com(str(self._id)+'USE'+str(num))
        self.com('USE'+str(num))
        #self.com('G')

    def del_profile(self,index=1):
        """ delete profile

        Parameters
        ----------

        index : int

        Examples
        --------

        >>> A.del_profile(num_profil)

        """
        self.lprofile.pop(index-1)

        assert(num<=len(self.lprofile)),"profile number not defined"
        self.com(str(self._id)+'USE'+str(num))


    def reset(self):
        """ reset axis

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.reset() #reset axis X

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

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes()
        >>> A.mvpro(1)

        """

        com = 'USE('+str(id_pro)+')'
        self.com(com)
        self.com('G')

    def read(self):
        pass

    def reg(self,typ='ST'):
        """ read boolean quantities in registers  : ST, UF, DF

        ST : Status  (default)
        UF : User Faults
        DF : Drive Faults

        Examples
        --------

        >>> # To check and parse the axis status
        >>> from pylayers.measures.parker import smparker
        >>> port = getty()
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.reg('ST')
        >>> #scans over axis 1 by given status

        """

        #buf = self.com('R','('+typ+')')
        buf = self.com('R'+'('+typ+')')
        buf = buf[1]
        buf = buf.replace('*','').replace('\r\n','').split('_')

        st =  ''
        for k in range(8):
            for l in range(4):
                val = eval(buf[k][l])
                if typ=='ST':
                    self.status[k*4+l] = val
                    if val:
                        st = st + ' ' + Axes.dstatus[k*4+l+1]+'\n'
                if typ=='UF':
                    self.usrflt[k*4+l] = val
                    if val:
                        st = st + ' ' + Axes.dusrflt[k*4+l+1]+'\n'
                if typ=='DF':
                    self.drvflt[k*4+l] = val
                    if val:
                        st = st + ' ' + Axes.ddrvflt[k*4+l+1]+'\n'
        return(st)

    def mv(self,var=0):
        """ move axes in translation or rotation

        Parameters
        ----------

        var : distance (cm) | degres (°)

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.mv(10) # moves over 10cm on axis 1
        >>> A.mv(45) # moves over 45° on axis 3

        """
        #assert(self.typ=='t'),'Axes is not a linear axes'
        #nstep = dcm*self.scale
        #com = self.com('D'+str(nstep))
        #if typ=='t':
            #nstep = dcm*self.scale
            #com = self.com('D'+str(nstep))

        nstep = int(var*self.scale)
        scom1 = 'D'+str(nstep)
        com = self.com(scom1)
        com = self.com('G')
        #com = self.com(scom1,verbose=True)
        #scom2 = 'G'
        #com = self.com(scom2,verbose=True)
        #print "distance parcourue : ", var+str('cm')
        #com = self.com(scom1,verbose=True)
        #com = self.com(scom2,verbose=True)
        #com = self.com(scom3,verbose=True)

    def close(self):
        self.ser.close()

    def util(self):
        """ allows convertion between :
            m/s | tr/s  <=> rps
            m/s²        <=> rps²
        """
        pass



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
        #st += self.ser.read(10)
        #self.ser.close()
        return(st)

class Scanner(PyLayers):
    """This class handles scenarios

    """

    def __init__(self,port=gettty()):
        """
        Parameters
        ----------

        p   : current position of the scanner
        phi : current angle of the scanner

        Examples
        --------

        >>> s.a[1].name_of_function()
        """
        self.ser = Serial(port = port, baudrate=9600, timeout = 1)
        #
        # p current position of the scanner
        #
        self.p = np.array([0,0])
        #
        # phi current angle of the scanner
        #
        self.phi = 0
        self.sx = 12800
        self.sy = 22800
        self.sr = 2111.111111111111
        self.a  = ['',Axes(1,'x',self.ser,scale=12800,typ='t'),
                      Axes(2,'y',self.ser,scale=22800,typ='t'),
                      Axes(3,'rot',self.ser,scale=2111.1111111111113,typ='r')]
                      #self.a4  = Axes(4,'z',self.ser,typ='r')
        self.a[1].limits(mode=1,typ=1,mask=0)
        self.a[2].limits(mode=1,typ=1,mask=0)
        self.a[3].limits(mode=1,typ=1,mask=0)

        self.a[1].com('ON')
        self.a[2].com('ON')
        self.a[3].com('ON')


    def __repr__(self):
        px = self.a[1].com('R(PA)')[1].replace('*','').replace('\n','')
        py = self.a[2].com('R(PA)')[1].replace('*','').replace('\n','')
        pr = self.a[3].com('R(PA)')[1].replace('*','').replace('\n','')
        st = str(float(px)/self.sx) +',' + str(float(py)/self.sy)  +','+ str(float(pr)/self.sr) +'\n'
        st = st + 'current position : '+ str(self.p) + '\n'
        st = st + 'current angle  : '+ str(self.phi) + '\n'
        return(st)


    def set_origin():
        """
        """
        pass

    def home(self,cmd='get'):
        """ allows a return home for 3 axes
        """
        for k in range(1,len(self.a)):
            self.a[k].home(cmd=cmd)

    def mv(pt,at,var=0):
        """ move to target point

        Parameters
        ----------

        pt : target position  (pt=np.array([0,0,0]))
        at : target angle

        """

        #
        #warphi : values prohibited to phi
        #
        warphi=np.arange(180,360)
        for i in warphi:
            assert(self.phi==i),'Error : Out of range'


        #Answers those questions
        # Ou suis-je ?
        # Ou dois-je aller : p1
        # Comment y aller :
        #   + fabriquer les profils
        #   + Appliquer les profils

    def array(A):
        """ Implement an Array




        """
        pass

        #array = ensemble de points
        #This fonction contains a cloud of points + noton of scheduling


if __name__=="__main__":
    port = gettty()
    X = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
    Y = Axes(2,'y',typ='t',scale=22800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
    X.limits(cmd='set',mode=1)
    Y.limits(cmd='set',mode=1)
    R = Axes(3,'r',typ='r',scale=2111.111111111111,ser=Serial(port=port,baudrate=9600,timeout=0.05))
    X.com()
    # pass
    #s = Scanner('/dev/ttyUSB0')
    #s = Scanner('/dev/ttyUSB2')
    #s = Scanner('/dev/ttyUSB1')
    #sm.fromfile('prog1')
    #sm.fromfile('AY')
    #Sc[1].com('ON')
    #st = sm.com(1,'LIMITS',(0,1,1))
    #st = sm.com(1,'1D-4000')
    #st = sm.com(1,'G')
