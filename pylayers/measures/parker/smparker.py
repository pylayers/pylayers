#!/usr/bin/python
#-*- coding:Utf-8 -*-
from serial import Serial
import pdb
import time
import threading
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pylayers.util.project import *
from pylayers.antprop.aarray import *
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

    if len(tty)>1:
        num = tty[1]
        port = '/dev/ttyUSB'+num
    else:
        port = None
        print 'not connected to a serial port'
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
        N       : number of point for plotting profile


        """

        defaults = {'num': 1,
                    'aa': 1,
                    'ad': 1,
                    'dstep': 16000,
                    'vmax': 5,
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
        self.tri = 2*np.sqrt(self.drev/self.aa)
        self.T = (  self.drev
                   + 0.5*self.vmax**2/self.aa
                   + 0.5*self.vmax**2/self.ad)/(1.0*self.vmax)
        if self.tri<self.T:
            self.T=self.tri

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
        self.t1 = min(self.vmax/(1.0*self.aa),self.T/2.0)
        self.t2 = max(self.T-self.vmax/(1.0*self.ad),self.T/2.0)

        u1 = np.where(self.t<self.t1)[0]
        u2 = np.where(self.t>=self.t2)[0]

        self.v[u1] = self.t[u1]*self.aa
        #self.v[u2] = -self.t[u2]*self.ad+(self.vmax+self.ad*t2)
        self.v[u2] = -self.t[u2]*self.ad+self.ad*self.T
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
        st = st + 'step per revolution  : '+ str(self.spr) + '\n'
        st = st + 'N  : '+ str(self.N) + '\n'
        return(st)

    def duration(self):
        """
        """
        #self.duration
        pass

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
        #for k in enumerate(str(self._id)):
            #st = st + '---------------------\n'
            #st = st + 'Info about motion '+' '+'Axe '+str(k[0]+1)+ '\n'
            #st = st + '---------------------\n'
            #st1 = self.reg()

        return(st)

    def show(self):
        """
        """
        pass

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
            st = self.com('R('+var+')')[1].replace('*','').replace('/n','')
            
            if var == 'MV':
                if '1' in st[1]:
                    print "I'm moving"
                else:
                    print "I'm stationnary"

            print Axes.dvar[var],st


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
            self.com(cstr,verbose=False)

    def home(self,cmd='get',**kwargs):
        """ enables back material home for each axe.
            for more informations see Parker book page 108.

        Parameters
        ----------

        cmd  : 'get','set','go'

        Examples
        --------


        """

        defaults = {'mode' :0,
                    'vel'  :15,
                    'acc'  :20,
                    'edg'  :'+',
                    'typ'  :0,
                    'armed':1
                }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.mode  = kwargs['mode']
        self.vel   = kwargs['vel']
        self.acc   = kwargs['acc']
        self.edg   = kwargs['edg']
        self.typ   = kwargs['typ']
        self.armed = kwargs['armed']

        if self._id in [1,2]:
            if cmd=='get':
                st = self.com('HOME')
                ans = st[1].split(' ')


                print "-------------------"
                print " arguments of HOME "
                print "-------------------"

                if '1' in ans[0]:
                    print "armed "
                else:
                    print "not armed "
                if '-' in ans[1]:
                    print "reference edge is negative "
                else:
                    print "reference edge is positive "
                if '1' in ans[2]:
                    print "home switch normally closed 1 "
                else:
                    print "home switch normally open 0 (default) "

                if '+' in ans[3]:
                    print 'velocity : +',eval(ans[3].split('V+')[1]), "rps"
                else:
                    print 'velocity : -',eval(ans[3].split('V-')[1]), "rps"

                print 'acceleration : ',eval(ans[4].split('A')[1]), "rps²"

                if '0' in ans[5]:
                    print 'Mode 0: Motor in the active window of the switch (default)'
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
                print cstr
                self.com(cstr)

            if cmd=='go':
                #cstr = 'HOME1'
                cstr = 'HOME'+str(self.armed)
                self.com(cstr)
                cstr = 'ARM'+str(self.armed)
                self.com(cstr)
                cstr = 'GH'
                #cstr = 'G'
                self.com(cstr)

        # no material origin available on z axisi(4) and rotation axis (3)
        if self._id in [3,4]:
            if cmd=='set':
                self.com('W(PA,0)')
            if cmd=='get':
                st =  self.com('R(PA)')
            if cmd=='go':
                pa = -int(self.com('R(PA)')[1].replace('*','').replace('\n',''))
                self.mv(pa/self.scale,vel=self.vel,acc=self.acc)

    def parsinfo(self,cmd='HOME'):
        """ allows get parsing informations. 

        Parameters
        ----------

        cmd  : 'LSEL' | 'HOME' | 'LIMITS'

        Examples
        --------


        """

        if cmd =='LSEL':
            st = self.com('LSEL')
            ans = st[1].split(' ')

            print "-------------------"
            print " arguments of LSEL "
            print "-------------------"

            ##########
            #### ARM 
            ##########
            if '1' in ans[0]:
                print "armed "
            else:
                print "not armed "

            ############ 
            ##### CODE
            ############

            if '1' in ans[1]:
                print "Binary (default setting) "
            else:
                print "BCD code (Binary Coded Decimal) "

            ############### 
            ##### INPUTS 
            ################

            if '5' in ans[2]:
                print "5 inputs "
            if '4' in ans[2]:
                print "4 inputs "
            if '3' in ans[2]:
                print "3 inputs "
            if '2' in ans[2]:
                print "2 inputs "
            if '1' in ans[2]:
                print "1 inputs "

            ##############
            #### EXECUTION 
            ##############

            if '0' in ans[3]:
                print "continuously repeated (default set.) "
            else:
                print "re-triggered "

        if cmd=='HOME':
            st = self.com('HOME')
            ans = st[1].split(' ')

            print "-------------------"
            print " arguments of HOME "
            print "-------------------"

            if '1' in ans[0]:
                print "armed "
            else:
                print "not armed "
            if '-' in ans[1]:
                print "reference edge is negative "
            else:
                print "reference edge is positive "
            if '1' in ans[2]:
                print "home switch normally closed 1 "
            else:
                print "home switch normally open 0 (default) "

            if '+' in ans[3]:
                print 'velocity : +',eval(ans[3].split('V+')[1]), "rps"
            else:
                print 'velocity : -',eval(ans[3].split('V-')[1]), "rps"
                print 'acceleration : ',eval(ans[4].split('A')[1]), "rps²"

            if '0' in ans[5]:
                print 'Mode 0: Motor in the active window of the switch (default)'
            if '1' in ans[5]:
                print 'Mode 1: Motor in the position to the edge + or -'
            if '2' in ans[5]:
                print 'Mode 2: Improve homing repeatability'

        if cmd=='LIMITS':
            st = self.com('LIMITS')
            ans = st[1].split(' ')

            print "-------------------"
            print " arguments of LIMITS "
            print "-------------------"


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


    def stationnary(self):
        """ test bit in status ST

        Examples
        --------

        >>> # To check and parse the axis status
        >>> from pylayers.measures.parker import smparker
        >>> port = getty()
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.reg('ST')

        >>> from pylayers.measures.parker import smparker
        >>> S = smparker.Scanner()
        >>> S.a[1].stationnary()
        >>> #scans over axis 1 by given status

        """

        buf = self.com('R'+'(ST)')
        buf = buf[1]
        buf = buf.replace('*','').replace('\r\n','').split('_')
        status = eval('0b'+reduce(lambda x,y:  x+y,buf))
        stationnary_bit = 0b00000000000000000001000000000000
        if eval(bin(stationnary_bit & status))!=0:
            return True
        else:
            return False


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

    def mv(self,var=0,vel=20,acc=30):
        """ move axes in translation or rotation

        Parameters
        ----------

        var : distance (cm) | degres (°)
        vel : velocity (rps)
        aa  : acceleration (rps²)

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
        >>> A.mv(10) # moves over 10cm on axis 1
        >>> A.mv(45) # moves over 45° on axis 3

        """
        nstep = int(var*self.scale) #convert num per step
        scom1 = 'D'+str(nstep) #command
        scom2 = 'V'+str(vel)   #set velocity
        scom3 = 'AA'+str(acc)   #set acceleration of motion
        #
        #send commands
        #
        t1 = time.time()
        com   = self.com(scom1)
        t2 = time.time()
        com   = self.com(scom2)
        t3 = time.time()
        com   = self.com(scom3)
        #tic = time.time()
        t4 = time.time()
        com = self.com('G')
        t5  = time.time()
        print t2-t1
        print t3-t2
        print t4-t3
        print t5-t4
        #toc = time.time()
        #print  toc-tic


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
    """This class handles the FACS (Four Axes Channel Scanner)

    """

    def __init__(self,port=gettty(),anchors={}):
        """
        Parameters
        ----------

        pG  : current position of the scanner in global frame
        pH  : current position of the scanner in home frame
        pA  : current position of the scanner in array frame
        ang : current angle of the scanner
        anchors : dict of scanner anchor points in global frame

        Examples
        --------

        >>> s.a[1].name_of_function()
        """
        self.ser = Serial(port = port, baudrate=9600, timeout = 1)
        self.anchors = anchors
        #
        # phi current angle of the scanner
        #
        self.phi = 0
        self.sx = 12800
        self.sy = 22800
        self.sz = 21111.1111111111
        self.sr = 2111.111111111111

        self.a  = ['',Axes(1,'x',self.ser,scale=self.sx,typ='t'),
                      Axes(2,'y',self.ser,scale=self.sy,typ='t'),
                      Axes(3,'rot',self.ser,scale=self.sr,typ='r'),
                      Axes(4,'z',self.ser,scale=self.sz,typ='t')]

        # p current position of the scanner
        # ang angle of the rotation axe
        #
        # Coordinate of Home Scanner Point in global frame
        if self.anchors=={}:
            self.H = np.array([0,0,0.9])
        else:
            pass
            #beware TBD from anchors

        # Coordinate of Array Scanner Point in home frame
        self.A = np.array([0,0,0.15])
        self.upd_pos(np.array([0,0,0]))
        self.ang =  0.
        # Limits activated on axes X and Y    (mask =0 )
        # Limits desactivated on axes Z and R (mask =3 )
        self.a[1].limits(mask=0,typ=1,mode=1,cmd='set')
        self.a[2].limits(mask=0,typ=1,mode=1,cmd='set')
        self.a[3].limits(mask=3,typ=1,mode=1,cmd='set')
        self.a[4].limits(mask=3,typ=1,mode=1,cmd='set')

        # Power ON
        self.a[1].com('ON')
        self.a[2].com('ON')
        self.a[3].com('ON')
        self.a[4].com('ON')

        self.home(cmd='set')
        self.home(cmd='go',init=True)


    def __repr__(self):
        st = 'Home frame : ' + str(self.pH[0])+','\
                                    + str(self.pH[1])+',' \
                                    + str(self.pH[2])+'\n'
        st = st +  'Global frame : '+str(self.pG[0])+','\
                                    + str(self.pG[1])+',' \
                                    + str(self.pG[2])+'\n'
        st = st + 'Array frame : '+str(self.pA[0])+','\
                                    + str(self.pA[1])+',' \
                                    + str(self.pA[2])+'\n'
        st =  st + 'Current angle : '+str(self.ang)+'\n'
        return(st)



    def check_pa(self):
        px = self.a[1].com('R(PA)')[1].replace('*','').replace('\n','')
        py = self.a[2].com('R(PA)')[1].replace('*','').replace('\n','')
        pz = self.a[4].com('R(PA)')[1].replace('*','').replace('\n','')

        pr = self.a[3].com('R(PA)')[1].replace('*','').replace('\n','')
        st = 'current position : '+ str(float(px)/self.sx) +',' + str(float(py)/self.sy)  +','+ str(float(pz)/self.sz) +'\n'
        st = st + 'current angle  : '+ str(float(pr)/self.sr) + '\n'


    def origin(self,cmd='set',p0=np.array([0,0,0]),ang0=0):
        """
        """
        if cmd=='set':
            self.p0=p0
            self.ang0=ang0

    def home(self,cmd='get',init=False,vel=10):
        """ allows a return home for 3 axes

        Parameters
        ----------

        cmd
        init
        vel
        frame

        """
        lt = []
        for k in range(1,len(self.a)):
            if init:
                if k in [3,4]:
                    cmd ='set'
            lt.append(threading.Thread(name=str(k),target=self.a[k].home(cmd=cmd,vel=10)))
        for k,t in enumerate(lt):
            t.start()
        self.upd_pos(np.array([0,0,0]))

    def upd_pos(self,ptH):
        self.pH = ptH
        self.pA = self.pH - self.A
        self.pG = self.H + self.pH

    def mv(self,pt=np.array([0.1,0.1,0]),at=0,frame='A',vel=20):
        """ move to target point

        Parameters
        ----------

        pt : target position (pt=np.array([0,0,0]))
        at : target angle
        frame : {'H'|'A'|'G'}

        """

        # convert to home frame
        if frame=='A':
            ptH = pt + self.A
        if frame=='G':
            ptH = pt + self.H
        if frame=='H':
            ptH = pt

        vec = ptH-self.pH
        dx = vec[0]
        dy = vec[1]
        dz = vec[2]
        da = at - self.ang
        self.a[1].mv(dx,vel=vel)
        self.a[2].mv(dy,vel=vel)
        self.a[4].mv(dz,vel=vel)
        self.a[3].mv(da,vel=vel)

        # update new position

        self.upd_pos(ptH)

        #Answers those questions:
        # Ou dois-je aller : p1
        # Comment y aller :
        #   + fabriquer les profils
        #   + Appliquer les profils

    def array(self,A,vel=25):
        """ Implement an Array

        Parameters
        ----------

        A : Aarray


        """
        for k in range(A.p.shape[1]):
            self.mv(pt=A.p[:,k],vel=vel)
            print self.a[1].stationnary()
            while not self.a[1].stationnary():
                print " I am moving"
            print "I am stationnary  now. You can measure"
            # while in motion
            #    pass
            # data =vna.measure()
            # thread store(data)
    def mes1(self,A,vna,vel=25):
        pass



        #array = ensemble de points
        #This fonction contains a cloud of points + noton of scheduling


if __name__=="__main__":
    S = Scanner()
    #vna =E()

    #S.a[axe]


    #run smparker
    #S=smparker.Scanner()
    #from pylayers.antprop.aarrray import *
    #A=AntArray()
    #S.array(A)


#    port = gettty()
#    X = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
#    Y = Axes(2,'y',typ='t',scale=22800,ser=Serial(port=port,baudrate=9600,timeout=0.05))
#    X.limits(cmd='set',mask=0)
#    Y.limits(cmd='set',mask=0)
#    R = Axes(3,'r',typ='r',scale=2111.111111111111,ser=Serial(port=port,baudrate=9600,timeout=0.05))
#    Z = Axes(4,'z',typ='t',scale=2111.111111111111,ser=Serial(port=port,baudrate=9600,timeout=0.05))
#    Z.limits(cmd='set',mask=3)
#    R.limits(cmd='set',mask=3)
#    # pass
#    #s = Scanner('/dev/ttyUSB0')
#    #s = Scanner('/dev/ttyUSB2')
#    #s = Scanner('/dev/ttyUSB1')
#    #sm.fromfile('prog1')
#    #sm.fromfile('AY')
#    #Sc[1].com('ON')
#    #st = sm.com(1,'LIMITS',(0,1,1))
#    #st = sm.com(1,'1D-4000')
#    #st = sm.com(1,'G')
