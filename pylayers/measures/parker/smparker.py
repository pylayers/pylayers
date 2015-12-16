#!/usr/bin/python
#-*- coding:Utf-8 -*-

from serial import Serial
import pdb
import time
import doctest
import threading
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pylayers.util.project import *
from pylayers.antprop.aarray import *
from pylayers.measures.exploith5 import *
from pylayers.measures.vna.E5072A import *
import pylayers.measures.switch.ni_usb_6501 as sw

"""

.. currentmodule:: pylayers.measures.parker

This module handles scanner.

Authors : M.D.BALDE and B.UGUEN

Profile Class
=============

.. autosummary::
    :toctree: generated/

    Profile.__init__
    Profile.__repr__
    Profile.duration
    Profile.show

Axes Class
==========

.. autosummary::
    :toctree: generated/

    Axes.__init__
    Axes.__repr__
    Axes.info
    Axes.show
    Axes.getvar
    Axes.com
    Axes.limits
    Axes.home
    Axes.parsinfo
    Axes.add_profile
    Axes.set_profile
    Axes.del_profile
    Axes.reset
    Axes.mvpro
    Axes.stationnary
    Axes.reg
    Axes.step
    Axes.velocity
    Axes.acceleration
    Axes.go
    Axes.close
    Axes.fromfile

Scanner Class
=============

.. autosummary::
    :toctree: generated/

    Scanner.__init__
    Scanner.__repr__
    Scanner.check_pa
    Scanner.origin
    Scanner.reset
    Scanner.home
    Scanner.upd_pos
    Scanner.mv
    Scanner.meash5

"""

def gettty():
    """get tty and handles port conflicts

    Examples
    --------

    >>> import  os
    >>> from pylayers.measures.parker import smparker

    """
    import  os
    line = os.popen('dmesg | grep tty | tail -1').read() .replace('\n','')
    tty = line.split('ttyUSB')
    if len(tty)>1:
        num = tty[1]
        port = '/dev/ttyUSB'+num
        return port

        if len(tty)>1:
            num = tty[1]
            port = '/dev/ttyUSB'+num
        else:
            port = None
            print 'not connected to a serial port'
    else:
        port = "not connected"

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
                    'aa': 5,
                    'ad': 1,
                    'dstep': 40000,
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
                 ser = 'emulated',
                 scale=1280000,
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
        >>> from pylayers.measures.parker.smparker import *
        >>> X = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> Y = Axes(2,'y',typ='t',scale=2280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> R = Axes(3,'r',typ='r',scale=2111.111111111111,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
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
        if self.ser == 'emulated':
            self.emulated = True
        else:
            self.emulated = False
    
        #
        # list of profiles
        #
        self.scale = scale
        self.typ = typ
        self.lprofile=[]
        self.add_profile()
        if not self.emulated:
            self.step()
            self.velocity()
            self.acceleration()

    def info(self):
        """gives informations about system variables of PARKER

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> st = A.info()
        """
        for k in self.dvar:
            st = self.com('R('+k+')')
            print self.dvar[k],st[1].replace('\n','')

    def __repr__(self):
        """
        """

        st = ''
        st = st + '---------------------\n'
        st = st + ' Parameters ' + '\n'
        st = st + '---------------------\n'
        st = st + 'axes_id : ' + str(self._id) + '\n'
        st = st + 'scale : ' + str(self.scale) + '\n'
        st = st + 'typ : ' + str(self.typ) + '\n'

        st  = st + '---------------------\n'
        st1 = self.reg()
        st  =  st + st1
        st  = st + '---------------------\n'

        if self.typ=='t':
            st =  st + 'dist : '+str(self.dist)+' m\n'
            st =  st + 'vel  : '+str(self.vel)+' rps\n'
            st =  st + 'acc  : '+str(self.acc)+' rps2\n'

        if self.typ=='r':
            st =  st + 'ang (deg) : '+str(self.ang)+' \n'
            st =  st + 'vel       : '+str(self.vel)+' rps\n'
            st =  st + 'acc       : '+str(self.acc)+' rps2\n'
        
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

    def getvar(self,lvar=[]):
        """allows get state of variables

        Parameters
        ----------

        lvar : list of variables

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
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
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> st=A.com('ON') #Energized motor
        >>> st=A.com('LIMITS(1,0,0)') #Disable limit +,limits normally open, mode 0 stop motion and abort prog
        >>> st=A.com('V15') # Change velocity to 15
        >>> st=A.com('AA15') # Change acceleration to 15
        >>> st=A.com('AD15') # Change deceleration to 15


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

        cmd   :  'get'|'set'
        mask  : int
            0 : Enable limits (default setting)
            1 : Disable limit +
            2 : Disable limit -
            3 : Disable limit + and -
        typ   : boolean
            0 : Limits normally closed (default setting)
            1 : Limits normally open
        mode  : boolean
            0 : Stop motion and abort the prog when a limit is hit
            1 : Stop motion and continue the prog when a limit is hit
        LD    : sets the required deceleration rate after hitting a limit

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

            print 'deceleration : ',eval(ans[3].split('D')[1]), "rps"

        if not self.emulated:
            if cmd=='set':
                cstr = 'LIMITS'+'('+str(mask)+','+str(typ)+','+str(mode)+','+str(LD)+')'
                self.com(cstr,verbose=False)

    def home(self,cmd='set',**kwargs):
        """ enables back material home for each axe.
            for more informations see Parker book page 108.

        Parameters
        ----------

        cmd   : 'set'|'go'
        mode  : int
            0 : The indexer positions the motor in the active window of the
                switch (default setting)
            1 : The motor is positioned to the required edge of the switch + or -
            2 : The indexer is set to seek the nearest drive zero phase
                position to improve homing repetability
        vel   : velocity and direction
        acc   : acceleration (and deceleration)
        edg   : this parameter is used to to select the required edge of the
                home switch
        typ   : boolean used to select the type of switch to be used for homing
            1 : Home switch normally closed
            0 : Home switch normally open (default)
        armed : boolean
            0 : turned OFF
            1 : turned ON

        Examples
        --------

        """

        defaults = {'mode' :0,
                    'vel'  :15,
                    'acc'  :25,
                    'edg'  :'+',
                    'typ'  :0,
                    'armed':1
                }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        mode  = kwargs['mode']
        vel   = kwargs['vel']
        acc   = kwargs['acc']
        edg   = kwargs['edg']
        typ   = kwargs['typ']
        armed = kwargs['armed']

        if self._id in [1,2]:
            if cmd=='set':
                if vel>0:
                    vel = '+'+str(vel)
                else:
                    vel = '-'+str(vel)

                cstr = 'HOME'+str(armed)+\
                              '('+edg+','+\
                              str(typ)+','+\
                              vel+','+\
                              str(acc)+','+\
                              str(mode)+')'
                self.com(cstr)

            if cmd=='go':
                cstr = 'HOME'+str(armed)
                self.com(cstr)
                cstr = 'ARM'+str(armed)
                self.com(cstr)
                cstr = 'GH'
                #cstr = 'G'
                self.com(cstr)

        # no material origin available on z axis(4) and rotation axis (3)
        if self._id in [3,4]:
            if cmd=='set':
                self.com('W(PA,0)')
            if cmd=='go':
                # to come back to origin the offset is -PA converted in the
                # proper scale
                pa = -int(self.com('R(PA)')[1].replace('*','').replace('\n',''))
                self.step(value=pa/(1.0*self.scale),cmd='set')
                self.go()

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

            ##############
            #### ARM
            ##############
            if '1' in ans[0]:
                print "armed "
            else:
                print "not armed "


            ##################
            ####REFERENCE EDGE
            ##################

            if '-' in ans[1]:
                print "reference edge is negative "
            else:
                print "reference edge is positive "



            ##############
            #### TYPE
            ##############

            if '1' in ans[2]:
                print "home switch normally closed 1 "
            else:
                print "home switch normally open 0 (default) "


            ##############################
            #### VELOCITY AND ACCELERATION
            ##############################
            if '+' in ans[3]:
                print 'velocity : +',eval(ans[3].split('V+')[1]), "rps"
            else:
                print 'velocity : -',eval(ans[3].split('V-')[1]), "rps"
                print 'acceleration : ',eval(ans[4].split('A')[1]), "rps"



            ##############
            #### MODE
            ##############

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


            ##############
            #### MASK
            ##############


            if '0' in ans[0]:
                print "Enable limits (default setting), "
            if '1' in ans[0]:
                print "Disable limit +, "
            if '2' in ans[0]:
                print "Disable limit -, "
            if '3' in ans[0]:
                print "Disable limit + & -, "

            ##############
            #### TYPE
            ##############

            if '0' in ans[1]:
                print "Limits normally closed (default setting), "
            else:
                print  "Limits normally open, "

            ##########################
            #### MODE AND DECELERATION
            ##########################

            if '0' in ans[2]:
                print "Stop motion when a limit is hit and abort the program (default setting), "
            else:
                print "Stop motion when a limit is hit but continue the program, "

            print 'deceleration : ',eval(ans[3].split('D')[1]), "rps"

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
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
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

    # def read(self):
    #     pass


    def stationnary(self):
        """ test bit in status ST

        Examples
        --------

        >>> # To check and parse the axis status
        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
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
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
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

    def step(self,value=0.1,cmd='get'):
        """  set/get distance

        Parameters
        ----------

        value : int
        cmd : string
            {set | get}

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> S = Scanner()
        >>> S.a[4].step(-0.1,cmd='set')
       

        """
        nstep = int(value*self.scale) #convert num per step
        if cmd=='set':
            scom = 'D'+str(nstep)   #set velocity
            if self.typ=='t':
                self.dist =  value
            else:
                self.ang = value
            if not self.emulated: 
                com   = self.com(scom)
        if cmd=='get':
            scom = 'D'   #get velocity
            rep =  self.com(scom)[1].replace('*','')
            if self.typ=='t':
                self.dist  = eval(rep)/(1.0*self.scale)
            else:
                self.ang  = eval(rep)/(1.0*self.scale)


    def velocity(self,value=10,cmd='get'):
        """  set/get velocity

        Parameters
        ----------

        value : int
        cmd : string
            {set | get}

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> A.velocity(value=10)
        >>> A.velocity(cmd='get')

        """
        if cmd=='set':
            scom = 'V'+str(value)   #set velocity
            self.vel =  value
            if not self.emulated:
                com   = self.com(scom)
        if cmd=='get':
            scom = 'V'   #get velocity
            self.vel  = eval(self.com(scom)[1].replace('*',''))


    def acceleration(self,value=10,cmd='get'):
        """  set/get acceleration

        Parameters
        ----------

        value : int
        cmd : string
            {set | get}

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> A.velocity(value=10)
        >>> A.velocity(cmd='get')

        """
        if cmd=='set':
            scom = 'AA'+str(value)   #set velocity
            self.acc =  value
            if not self.emulated:
                com   = self.com(scom)
        if cmd=='get':
            scom = 'AA'   #get velocity
            self.acc  = eval(self.com(scom)[1].replace('*',''))


    def go(self):
        """ move axes in translation or rotation

        Parameters
        ----------

        var : distance (cm) | degres ()
        vel : velocity (rps)
        aa  : acceleration (rps)

        Examples
        --------

        >>> from pylayers.measures.parker import smparker
        >>> A = Axes(1,'x',typ='t',scale=1280000,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> R = Axes(3,'ang',typ='r',scale=2111.1111,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
        >>> A.step(value=0.1,cmd='set')
        >>> A.go() # moves over 10cm on axis 1
        >>> R.go() # moves over 45 on axis 3

        """
        if not self.emulated:
            com = self.com('G')
            while not self.stationnary():
                pass


    def close(self):
        self.ser.close()

    # def util(self):
    #     """ allows convertion between :
    #         m/s | tr/s  <=> rps
    #         m/s        <=> rps
    #     """
    #     pass



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
    """ This class handles the FACS (Four Axes Channel Scanner)

    """

    def __init__(self,port=gettty(),anchors={},reset=True,vel=15,acc=15,**kwargs):
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

        defaults = {'time':True,
                    'Nt' : 4,
                    'Nr' : 8,
                    'emulated' :True}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        time = kwargs.pop('time')
        self.emulated   = kwargs.pop('emulated')
        self.Nt   = kwargs.pop('Nt')
        self.Nr   = kwargs.pop('Nr')

        if not self.emulated:
            self.ser = Serial(port = port, baudrate=9600, timeout = 1)
        else:
            self.ser = 'emulated'
        self.anchors = anchors


        #
        # phi current angle of the scanner
        #
        self.phi = 0
        #
        #  scale is expressed in step/m for translation axis 1,2,4
        #  and in step/deg for rotation axis 3
        #
        #  alpha : m/tour
        #  beta  : reduction
        #  N     : step/tour
        #  N/alpha  : step/m
        #,
        #  Axis 1 : 0.003125 m/tr 1   4000 step/tr 1280000 step/m
        #  Axis 2 : 0.005 m/tr  2.85  4000 step/tr 2280000 step/m
        #  Axis 4 : 0.004 m/tr  1     4000 step/tr 1000000 step/m
        #
        self.sx = 1280000
        self.sy = 2280000
        self.sz = 1000000
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
            self.H = np.array([0,0,0])
        else:
            self.H = (anchors['p1']+anchors['p2']+anchors['p3'])/3.
            #beware TBD from anchors

        # Coordinate of Array Scanner Point in home frame
        self.A = np.array([0,0,0.1])
        self.upd_pos(np.array([0,0,0]))
        self.ang =  0.
        # Limits activated on axes X and Y    (mask =0 )
        # Limits desactivated on axes Z and R (mask =3 )

        #if not self.emulated:
        print "setting limits"
        self.a[1].limits(mask=0,typ=1,mode=1,cmd='set')
        self.a[2].limits(mask=0,typ=1,mode=1,cmd='set')
        self.a[3].limits(mask=3,typ=1,mode=1,cmd='set')
        self.a[4].limits(mask=3,typ=1,mode=1,cmd='set')

    #print "reseting axes"
    #if reset:
    #    self.reset()


        print "setting step"
        self.a[1].step(0,cmd='set')
        self.a[2].step(0,cmd='set')
        self.a[3].step(0,cmd='set')
        self.a[4].step(0,cmd='set')

        print "setting velocity"
        self.a[1].velocity(vel,cmd='set')
        self.a[2].velocity(vel,cmd='set')
        self.a[3].velocity(vel,cmd='set')
        self.a[4].velocity(vel,cmd='set')

        print "setting acceleration"
        self.a[1].acceleration(acc,cmd='set')
        self.a[2].acceleration(acc,cmd='set')
        self.a[3].acceleration(acc,cmd='set')
        self.a[4].acceleration(acc,cmd='set')

        #self.home(cmd='set',vel=10) #vel = 1 by default
        #print "home"
        #self.home(cmd='go',init=True)


    def __repr__(self):
        """
        """

        st = ''
        st = st + '---------------------------\n'
        st = st + '    Parameters of Scan' + '\n'
        st = st + '---------------------------\n'
        st = st + 'Home frame   : ' + str(self.pH[0])+','\
                                    + str(self.pH[1])+',' \
                                    + str(self.pH[2])+'\n'
        st = st + 'Global frame : '+str(self.pG[0])+','\
                                    + str(self.pG[1])+',' \
                                    + str(self.pG[2])+'\n'
        st = st + 'Array frame  : '+str(self.pA[0])+','\
                                    + str(self.pA[1])+',' \
                                    + str(self.pA[2])+'\n'
        st = st + 'Current angle : '+str(self.ang)+'\n'

        if not self.emulated:
            st = st + '-------------------\n'
            st = st + 'x (d,v,a)     :'+ str(self.a[1].dist)+' '+str(self.a[1].vel)+' '+str(self.a[1].acc)+'\n'
            st = st + 'y (d,v,a)     :'+ str(self.a[2].dist)+' '+str(self.a[2].vel)+' '+str(self.a[2].acc)+'\n'
            st = st + 'z (d,v,a)     :'+ str(self.a[4].dist)+' '+str(self.a[4].vel)+' '+str(self.a[4].acc)+'\n'
            st = st + '------------------\n'
            st = st + 'rot (a,w,w2)  : '+ str(self.a[3].ang)+' ' +str(self.a[3].vel)+' '+str(self.a[3].acc)+'\n'

        return(st)



    def check_pa(self):
        """
        """
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

    def reset(self):
        """ reset and enpower all axes
        """
        tic = time.time()
        for k in range(1,len(self.a)):
            com = self.a[k].reset()
        toc = time.time()
        print "time reset (s) :",toc-tic


    def home(self,cmd='set',init=True,vel=10):
        """ allows a return home for 3 axes

        Parameters
        ----------

        cmd   : set
        init  : boolean (False)
        vel   : velocity (10)
        frame : landmark {'H'|'A'|'G'}

        """
        for k in range(1,len(self.a)):
            if init:
                if k in [3,4]:
                    cmd ='set'
            self.a[k].home(cmd=cmd,vel=10)
        self.upd_pos(np.array([0,0,0]))

    def upd_pos(self,pH):
        """ update position

        Parameters
        ----------

        pH : current position in H frame

        Returns
        -------

        self.pH : corrdinate in Home frame
        self.pA : coordinate in Array frame
        self.pG : coordinate in Global Frame

        """

        self.pH = pH
        self.pA = self.pH - self.A
        self.pG = self.H + self.pH

    def mv(self,pt=np.array([0.,0.,0]),at=0,frame='A',vel=20):
        """ move to target point

        Parameters
        ----------

        pt : target position (pt=np.array([0,0,0]))
        at : target angle
        frame : {'G'|'H'|'A'}
            determine in which frame is expressed pt (default A)
        vel : int
            velocity 
        """

        # convert to home frame
        # pt : target po
        # ptH : target point in Home Frame
        # self.A : [0,0,0.1]_H it means that the array origin is 10 cm above Home frame origin
        # self.H : [10,10,1.5]_G is the position of origin of Home Frame scanner in Global frame
        #
        if frame == 'A':  # pt expressed in A
            ptH = pt + self.A
        if frame == 'G':  # pt expressed in G
            ptH = pt + self.H
        if frame == 'H':
            ptH = p0
            #ptH = pt

        # self.pH : current position in Home frame
        # ptH     : target point in Home Frame

        vec = ptH - self.pH
        dx = vec[0]
        dy = vec[1]
        dz = vec[2]
        da = at - self.ang
        #print "source : ",self.pH
        #print "target : ",ptH
        #print "dx,dy,dz,da : ",dx,dy,dz,da
        #print "a1,a2,a4,a3:",self.a[1].dist,self.a[2].dist,self.a[4].dist,self.a[3].ang


        # move axis only if there is a modification from the previous move
        if dx!=0:
            if dx!=self.a[1].dist:
                self.a[1].step(value=dx,cmd='set')
            self.a[1].go()
        if dy!=0:
            if dy!=self.a[2].dist:
                self.a[2].step(value=dy,cmd='set')
            self.a[2].go()
        if dz!=0:
            if dz!=self.a[4].dist:
                self.a[4].step(value=dz,cmd='set')
            self.a[4].go()
        if da!=0:
            if da!=self.a[4].ang:
                self.a[3].step(value=da,cmd='set')
            self.a[3].go()

        # update new position
        self.upd_pos(ptH)


    def meash5(self,
               A,
               _fileh5='test.h5',
               gcal = 1,
               ical = 1,
               vel = 15,
               Nmeas = 1,
               pAnt = np.array([1.6,5.2,1.6]),
               comment = '',
               author = ''):
        """ measure over a set of points from AntArray and store in h5

        Parameters
        ----------

        A       : Aarray
        _fileh5 : string
            name of the h5 file containing calibration data
        gcal    : calibration group
        ical    : calibration index
        vel     : int
            scanner moving velocity
        Nmeas   : int
            Number of measurement
        Nr : int
            Number of receiving antenna
        Nt : int
            Number of transmitting antenna
        pAnt : np.array(,3)
            Coordinates of antenna phase center

        Examples
        --------

        """
        # load the file containing the calibration data
        if '.h5' not in _fileh5:
            _fileh5 = _fileh5+'.h5'

        Dh5 = mesh5(_fileh5)
        # open - sdata analysis
        Dh5.open('r')
        try:
            ldataset = Dh5.f.keys()
        except:
            raise IOError('no calibration in h5 file')
        lcal= [ eval(k.replace('cal','')) for k in ldataset if 'cal' in k ]
        lcal=np.array(lcal)

        if len(lcal)==1:
           gcal = lcal[0]
        else:
            if gcal not  in lcal:
                raise IOError('Error calibration : File does not exist')
        Dh5.close()
        # read the chosen calibration and save parameters in ini file for VNA

        Dh5.readcal(gcal=gcal,ical=ical)
        # update vna_config.ini
        Dh5.saveini()
        # end of read and save
        # initialization of vna
        vna = SCPI()
        vna.load_config_vna()

        Npoint = A.p.shape[1]
        Nf = vna.Nf

        laxes = []
        if A.N[0]!=1:
            laxes.append('x')
        if A.N[1]!=1:
            laxes.append('y')
        if A.N[2]!=1:
            laxes.append('z')
        if A.N[3]!=1:
            laxes.append('a')
        lN =  [ A.N[k] for  k  in range(3) if A.N[k]!=1 ]

        #
        # # end of initialization

        Dh5.open('a')
        try:
            ldataset = Dh5.f.keys()
        except:
            ldataset = []
        lmes = [ldataset[k] for  k in range(len(ldataset))  if 'mes' in ldataset[k]]
        mesname = 'mes'+str(len(lmes)+1)

        mes = Dh5.f.create_group(mesname)

        mes.attrs['time'] = time.ctime()
        mes.attrs['author'] = author
        mes.attrs['comment'] = comment
        mes.attrs['axes'] = laxes
        mes.attrs['axesN'] = lN
        mes.attrs['Nr'] = self.Nr
        mes.attrs['Nt'] = self.Nt
        mes.attrs['Nmeas'] = Nmeas
        mes.attrs['pAnt'] = pAnt
        #mes.attrs['anchors']= self.anchors
        # here is the hard link between a measurement and its calibration 
        mes.attrs['gcal'] = "cal"+str(gcal)
        mes.attrs['ical'] = str(ical)


        # Measure
        ik = np.arange(A.p.shape[1])
        ix,iy,iz,ia = k2xyza(ik,(A.N[0],A.N[1],A.N[2],A.N[3]))

        # k iterates on the total number of points

        for k in ik:
            self.mv(pt=A.p[:,k],vel=vel)
            # call vna for measurement

            S = vna.getdata(Nmeas=Nmeas,Nr=self.Nr,Nt=self.Nt)

            # if  len(lN)==1:
            #     mes.create_dataset(str(k),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
            #     mes[str(k)].attrs['pt']=A.p[:,k]
            #     mes[str(k)].attrs['pG']=self.pG
            #     mes[str(k)].attrs['pA']=self.pA

            # if  len(lN)==2:
            #     i2,i1,i0 = ktoxyz(k,N1=lN[1],N2=lN[0])
            #     try:
            #         mes.create_group(str(i1+1))
            #         mes[str(i1)].create_dataset(str(i0),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
            #     except:
            #         mes[str(i1)].create_dataset(str(i0),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
            #     mes[str(i1)][str(i0)].attrs['pt']=A.p[:,k]
            #     mes[str(i1)][str(i0)].attrs['pG']=self.pG
            #     mes[str(i1)][str(i0)].attrs['pA']=self.pA

            # if  len(lN)==3:
            try:
                mes.create_group(str(ix[k]))
                try:
                    mes[str(ix[k])].create_group(str(iy[k]))
                    try: 
                        mes[str(ix[k])][str(iy[k])].create_group(str(iz[k]))
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
                    except:
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
                except:
                    try: 
                        mes[str(ix[k])][str(iy[k])].create_group(str(iz[k]))
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
                    except:
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
            except:
                try:
                    mes[str(ix[k])].create_group(str(iy[k]))
                    try: 
                        mes[str(ix[k])][str(iy[k])].create_group(str(iz[k]))
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
                    except:
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
                except:
                    try: 
                        mes[str(ix[k])][str(iy[k])].create_group(str(iz[k]))
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)
                    except:
                        mes[str(ix[k])][str(iy[k])][str(iz[k])].create_dataset(str(ia[k]),(Nmeas,self.Nr,self.Nt,Nf),dtype=np.complex64, data = S)

            mes[str(ix[k])][str(iy[k])][str(iz[k])][str(ia[k])].attrs['pt']=A.p[:,k]
            mes[str(ix[k])][str(iy[k])][str(iz[k])][str(ia[k])].attrs['pG']=self.pG  
            mes[str(ix[k])][str(iy[k])][str(iz[k])][str(ia[k])].attrs['pA']=self.pA    

        Dh5.close()
#
#    def meash5(self,
#               A,
#               _fileh5='test.h5',
#               gcal=1,
#               ical=1,
#               vel=15,
#               Nmeas=1,
#               comment='',
#               author=''):
#        """ measure over a set of points from AntArray and store in h5
#
#        Parameters
#        ----------
#
#        A       : Aarray
#        _fileh5 : string
#            name of the h5 file containing calibration data
#        ical    : calibration group
#        vel     : int
#            scanner moving velocity
#        Nmeas   : int
#            Number of measurement
#
#        Examples
#        --------
#
#        """
#
#        # load the file containing the calibration data
#        if '.h5' not in _fileh5:
#            _fileh5 = _fileh5+'.h5'
#
#        Dh5 = mesh5(_fileh5)
#        # open - sdata analysis
#        Dh5.open('r')
#        try:
#            ldataset = Dh5.f.keys()
#        except:
#            raise IOError('no calibration in h5 file')
#        lcal= [ eval(k.replace('cal','')) for k in ldataset if 'cal' in k ]
#        lcal=np.array(lcal)
#
#        if len(lcal)==1:
#           ical = lcal[0]
#        else:
#            if ical not  in lcal:
#                raise IOError('Error calibration : File does not exist')
#        Dh5.close()
#        # read the chosen calibration and save parameters in ini file for VNA
#        Dh5.readcal(gcal=gcal,ical=ical)
#        # update vna_config.ini
#        Dh5.saveini()
#        # end of read and save
#        # initialization of vna
#        vna = SCPI()
#        vna.load_config_vna()
#
#        Npoint = A.p.shape[1]
#        laxes = []
#        if A.N[0]!=1:
#            laxes.append('x')
#        if A.N[1]!=1:
#            laxes.append('y')
#        if A.N[2]!=1:
#            laxes.append('z')
#        lN =  [ A.N[k] for  k  in range(3) if A.N[k]!=1 ]
#        
#        pdb.set_trace()
#        Nf = vna.Nf
#        # end of initialization
#
#        Dh5.open('a')
#        try:
#            ldataset = Dh5.f.keys()
#        except:
#            ldataset = []
#        lmes = [ldataset[k] for  k in range(len(ldataset))  if 'mes' in ldataset[k]]
#        mesname = 'mes'+str(len(lmes)+1)
#
#        if  laxes==['x']:
#            mes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],Nf),dtype=np.complex64)
#
#        if  laxes==['y']:
#            mes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],Nf),dtype=np.complex64)
#
#        if  laxes==['z']:
#            mes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],Nf),dtype=np.complex64)
#
#        if  laxes==['x','y']:
#            mes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],lN[1],Nf),dtype=np.complex64)
#
#        if  laxes==['x','z']:
#            mes =  Dh5.f.create_dataset(mesname,(Nmeas,lN[0],lN[2],Nf),dtype=np.complex64)
#
#        if  laxes==['y','z']:
#            mes = Dh5.f.create_dataset(mesname,(Nmeas,lN[1],lN[2],Nf),dtype=np.complex64)
#
#        if  laxes==['x','y','z']:
#            mes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],lN[1],lN[2],Nf),dtype=np.complex64)
#
#        mes.attrs['time'] = time.ctime()
#        mes.attrs['author'] = author
#        mes.attrs['comment'] = comment
#        mes.attrs['axes'] = laxes
#        #mes.attrs['min'] = lmin
#        #mes.attrs['max'] = lmax
#        mes.attrs['Nmeas'] = Nmeas
#        # here is the hard link between a measurement and its calibration 
#        mes.attrs['gcal'] = "cal"+str(gcal)
#        mes.attrs['ical'] = str(ical)
#
#
#        # Measure
#        #A.p.shape : Naxis x Npoint (3,8)
#
#        for k in np.arange(A.p.shape[1]):
#            # A.p[:,k].shape : (3,)
#            self.mv(pt=A.p[:,k],vel=vel)
#            S = vna.getdata(Nmeas=Nmeas)
#            if  laxes==['x']:
#                mes[:,k,:] = S
#
#            if  laxes==['y']:
#                mes[:,k,:] = S
#
#            if  laxes==['z']:
#                mes[:,k,:] = S
#
#            if  laxes==['x','y']:
#                ix,iy = ktoxyz(k,Nx=l)
#                mes[:,ix,iy,:] = S
#
#            if  laxes==['x','z']:
#                ix,iz = ktoxyz(k,Nx=l)
#                mes[:,ix,iz,:] = S
#
#            if  laxes==['y','z']:
#                iy,iz = ktoxyz(k,Nx=l)
#                mes[:,iy,iz,:] = S
#
#            if  laxes==['x','y','z']:
#                ix,iy,iz = ktoxyz(k,Nx=l)
#                mes[:,ix,iy,iz,:] = S
#
#        Dh5.close()
    def measMIMO(self,
               A,
               _fileh5='test.h5',
               ical=1,
               vel=15,
               Nmeas=100,
               comment='',
               author='',
               **kwargs):
        """ measure MIMO channel over a set of point from AntArray and store in h5

        Parameters
        ----------

        A       : Aarray
        _fileh5 : string
            name of the h5 file containing calibration data
        vel     : int
            scanner moving velocity
        Nmeas   : int
            Number of measurement

        """
  
        #initialization of the switch

        switch = sw.get_adapter()
        reattach=False
        if not switch:
            raise Exception("No device found")

        switch.device
        switch.set_io_mode(0b11111111, 0b11111111, 0b00000000)
        
        # load the file containing the calibration data
        Dh5 = mesh5(_fileh5)
                
        # open - sdata analysis
        Dh5.open('r')
        
        try:
            ldataset = Dh5.f.keys()
        except:
            raise IOError('no calibration in h5 file')

        lmimocal= [ eval(k.replace('cal','')) for k in ldataset if 'cal' in k ]
        lmimocal=np.array(lmimocal)

        if len(lmimocal)==1:
           imimocal = lmimocal[0]
        else:
            if imimocal not in lmimocal:
                raise IOError('Error calibration MIMO : File does not exist')
        Dh5.close()

        # read the chosen calibration and save parameters in ini file for VNA

        Dh5.readcal(imimocal=imimocal)
        Dh5.saveini()

        # end of read and save

        # initialization of vna

        vna = SCPI()
        vna.load_config_vna()

        Npoint = A.p.shape[1]
        laxes = []
        if A.N[0]!=1:
            laxes.append('x')
        if A.N[1]!=1:
            laxes.append('y')
        if A.N[2]!=1:
            laxes.append('z')

        lN =  [ A.N[k] for  k  in range(3) if A.N[k]!=1 ]
        Nf = vna.Nf

        # end of initialization

        Dh5.open('a')

        try:
            ldataset = Dh5.f.keys()
        except:
            ldataset = []

    
        # lmimomes = [ldataset[k] for  k in range(len(ldataset))  if 'mes' in ldataset[k]]
        # mesname = 'mes'+str(len(lmimomes)+1)

        # mimomes = Dh5.f.create_dataset(mesname,(Nmeas,self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['x']:
        #     mimomes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['y']:
        #     mimomes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['z']:
        #     mimomes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['x','y']:
        #     mimomes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],lN[1],self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['x','z']:
        #     mimomes =  Dh5.f.create_dataset(mesname,(Nmeas,lN[0],lN[2],self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['y','z']:
        #     mimomes = Dh5.f.create_dataset(mesname,(Nmeas,lN[1],lN[2],self.Nt,self.Nr,Nf),dtype=np.complex64)

        # if  laxes==['x','y','z']:
        #     mimomes = Dh5.f.create_dataset(mesname,(Nmeas,lN[0],lN[1],lN[2],self.Nt,self.Nr,Nf),dtype=np.complex64)


        # mimomes.attrs['time'] = time.ctime()
        # mimomes.attrs['author'] = author
        # mimomes.attrs['comment'] = comment
        # mimomes.attrs['axes'] = laxes
        # mimomes.attrs['Nmeas'] = Nmeas

        # # here is the hard link between a measurement and its calibration

        # for iR in range(self.Nr):
        #     for iT in range(self.Nt):
        #         mimomes.attrs['mimocal'] = "mimocal"+str(imimocal)

        # # Measure
        # #A.p.shape : Naxis x Npoint (3,8)

        # for k in np.arange(A.p.shape[1]):
        #     self.mv(pt=A.p[:,k],vel=vel)
        #     for iR in range(self.Nr):
        #         switch.write_port(0,iR)
        #         for iT in range(self.Nt):
        #             print iR,iT
        #             switch.write_port(1,iT)
        #             Smeas = vna.getdata(Nmeas=Nmeas)

        #             if  laxes==['x']:
        #                 mimomes[:,k,:] = Smeas

        #             if  laxes==['y']:
        #                 mimomes[:,k,:] = Smeas

        #             if  laxes==['z']:
        #                 mimomes[:,k,:] = Smeas

        #             if  laxes==['x','y']:
        #                 ix,iy = ktoxyz(k,Nx=l)
        #                 mimomes[:,ix,iy,:] = Smeas

        #             if  laxes==['x','z']:
        #                 ix,iz = ktoxyz(k,Nx=l)
        #                 mimomes[:,ix,iz,:] = Smeas

        #             if  laxes==['y','z']:
        #                 iy,iz = ktoxyz(k,Nx=l)
        #                 mimomes[:,iy,iz,:] = Smeas

        #             if  laxes==['x','y','z']:
        #                 ix,iy,iz = ktoxyz(k,Nx=l)
        #                 mimomes[:,ix,iy,iz,:] = Smeas

        Dh5.close()


    # def meas(self,A,vel=10,Nmeas=1):
    #     """ Measure over a set of point from AntArray

    #     Parameters
    #     ----------

    #     A : Aarray


    #     """
    #     vna = SCPI()
    #     vna.load_config()

    #     Npoint = A.p.shape[1]
    #     Nf = vna.Nf
    #     Smeas = np.empty((Nmeas,Npoint,Nf),dtype=complex)
    #     print Smeas.shape
    #     for k in np.arange(A.p.shape[1]):
    #         print k
    #         print A.p[:,k]
    #         # find a rule to retrieve ix,iy,iz,ia from k
    #         self.mv(pt=A.p[:,k],vel=vel)
    #         # Nmeas x Nf
    #         S = vna.getdata(Nmeas=Nmeas)
    #         Smeas[:,k,:]=S.y
    #     return(Smeas)

if __name__=="__main__":
    doctest.testmod()
    #S = Scanner()
    #vna =E()

    #S.a[axe]


    #run smparker
    #S=smparker.Scanner()
    #from pylayers.antprop.aarrray import *
    #A=AntArray()
    #S.array(A)


#    port = gettty()
#    X = Axes(1,'x',typ='t',scale=12800,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
#    Y = Axes(2,'y',typ='t',scale=22800,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
#    X.limits(cmd='set',mask=0)
#    Y.limits(cmd='set',mask=0)
#    R = Axes(3,'r',typ='r',scale=2111.111111111111,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
#    Z = Axes(4,'z',typ='t',scale=2111.111111111111,ser=Serial(port=gettty(),baudrate=9600,timeout=0.05))
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