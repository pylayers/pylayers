#!/usr/bin/python
#-*- coding:Utf-8 -*-
from serial import Serial
import pdb
import numpy as np
import matplotlib.pyplot as plt


class Profile(object):
    accmax = 200
    vmax = 10
    dmax = 1500

    def __init__(self,num,aa,ad,dstep,vtri,vs=0,spr=4000):
    #def __init__(self,num,aa,ad,dstep,v,vs=0,spr=4000):
        """

        Parameters
        ----------

        num     : profile number 
        aa      : aceleration (rps2)
        ad      : deceleration (rps2)
        dstep   : distance (step)
        v       : vmax (rps)
        vs      : vstart (rps)
        spr     : steps per round

        """ 

        N           = 100
        self.dstep  = dstep
        #self.dround = dstep/(1.0*spr)
        #self.T      = self.dround/(1.0*v)
        
        self.droundtri = dstep/(1.0*spr)
        self.Ttri      = self.droundtri/(1.0*vtri)
        
        #assert(0<aa<accmax)
        #assert(0<v<vmax)
        #assert(0<d<dmax)

        #self.cmd = 'PROFILE'+str(num)+'('+str(aa)+','+str(ad)+','+str(dstep)+','+str(v)+','+str(vs)+')'
        self.cmd = 'PROFILE'+str(num)+'('+str(aa)+','+str(ad)+','+str(dstep)+','+str(vtri)+','+str(vs)+')'
        
        #self.t     = np.linspace(0,self.T,N)
        #self.v     = v*np.ones(N)
        #self.v1    = self.t*aa
        #self.v2    = self.t*ad
        #t1         = v/(1.0*aa)
        #t2         = self.T-v/(1.0*ad)
        #u1         = np.where(self.t<t1)[0]
        #u2         = np.where(self.t>=t2)[0]
        #self.v[u1] = self.t[u1]*aa
        #self.v[u2] = -self.t[u2]*ad+(v+ad*t2)
        #self.dr    = np.cumsum(self.v)*(self.t[1]-self.t[0])
        #self.ds    = self.dr*spr
        
        #PROFILE MODELE TRIANGLE
        self.ttri   = np.linspace(0,self.Ttri,N)
        self.vtri   = vtri*np.ones(N)
        
        self.vtri1  = self.ttri*aa
        #self.vtri2  = self.ttri*ad
        
        ttri1       = vtri/(1.0*aa)
        #ttri2       = -vtri/(1.0*ad)

        utri1 = np.where(self.ttri<ttri1)[0]
        utri2 = np.where(self.ttri>=ttri1)[0]
        
        self.vtri[utri1]   = self.ttri[utri1]*aa
        #self.vtri[utri2]   = self.ttri[utri2]*ad
        self.vtri[utri2]   = self.ttri[utri1]*ad


        self.droundtri  = np.cumsum(self.vtri)*(self.ttri[1]-self.ttri[0])
        self.dstri      = self.droundtri*spr
        
        

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
        #plt.plot(self.t,self.v)
        plt.plot(self.ttri,self.vtri)
        plt.xlabel('time (s)')
        plt.ylabel('Velocity (rev/s)')
        plt.title('Evoltion of velocity over time')
        plt.subplot(212)
        #plt.plot(self.ds,self.v)
        plt.plot(self.dstri,self.vtri)
        plt.xlabel('distance(rev)')
        plt.ylabel('Velocity (rev/s)')
        plt.title('Evoltion of velocity over distance')
        plt.legend()
        plt.show()

    

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
    ddrvflt[3]='Supply rail failure'
    ddrvflt[4]='Ambient over temperature'
    ddrvflt[5]='Drive over temperature'
    ddrvflt[6]='Configuration error'
    ddrvflt[7]='Motor high voltage rail failure'
    ddrvflt[8]='Output fault'

    def __init__(self,_id,name,ser,scale=12800,typ='t'):
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
        # serial port
        self.ser = ser
        # list of profiles
        self.scale = scale
        self.typ = typ
        self.lprofile=[]



    def show(self):
        """
        """
        pass

    #def com2(self,com='R(ST)'):
        #self.ser.write(str(self._id)+com+'\r\n')
        #st = self.ser.readlines()
        #return(st)

    def com(self,name,rg='',verbose=False):
        if rg!='':
            cst = str(self._id)+name+str(rg)+'\r\n'
        else:
            cst = str(self._id)+name+'\r\n'
        if verbose:
            print cst
        self.ser.write(cst)
        st = self.ser.readlines()
        return(st)


   
    def home(self):
        """ enables back home
        """
        #if axis!=[]:
        if str(self._id)!=[]:
            cstr = 'HOME1(-,1,-15,100,2)'
        else:
            cstr = 'HOME1(-,1,-15,100,2)'
        self.com(cstr)
        self.com('GH')
        err = self.reg('UF')
        return(err)

    def add_profile(self,aa,ad,d,v,vs):
        """  Add new profile to list
        """
        npro = len(self.lprofile)
        if len(self.lprofile<8):
            prof = Profile(npro+1,aa,ad,d,v,vs)
            # update profile
            self.com(prof.cmd)
            self.lprofile.append(prof)
        
    def set_profile(self,num):
        """
        """
        assert(num<=len(self.lprofile)),"profile number not defined" 
        self.com(str(self._id)+'USE'+str(num))   

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

        sm.mvpro(1,1)  # axis 1 , profile 1
        """

        #com = 'USE('+str(id_pro)+')'
        #self.com(com)
        #self.com('G')
        
        com = 'USE('+str(id_pro)+')'
        self.com(com)
        self.com('G')

    def read(self):
        pass

    def reg(self,typ='ST'):
        """ read boolean quantities in registers  : ST,UF,DF
    
        Examples
        --------
          
        s.a[1].reg('ST') #scans over axis 1 by given status
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
                        print Axes.dusrflt[k*4+l+1]
                if typ=='DF':
                    self.drvflt[k*4+l] = val
                    if val:
                        print Axes.ddrvflt[k*4+l+1]

    def mv(self,var=0):
        """ move axes in translation or rotation
        
        Parameters
        ----------
        
        var : distance (cm) | degres (°)
       
        Examples
        --------

        s.a[1].mv(10)  # moves over 10cm on axis 1
        s.a[3].mv(45) # moves over 45° on axis 3
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
           
    def mvhome(self):
        """Enables going home (last position where it was)
        """
        #pass
        #nstep = int(var*self.scale)
        #scom1 = 'D'+str(nstep)
        scom2 = 'H'
        scom3 = 'G'
        
        #com = self.com(scom1,verbose=True)
        com = self.com(scom2,verbose=True)
        com = self.com(scom3,verbose=True)
        
    def homeor(self):
        """Back to material origin
        """
        
        scom0 = 'HOME'+ str(self._id)
        scom1 = 'ARM'+ str(self._id)
        com = self.com(scom0)
        com = self.com(scom1)
        com = self.com('GH')
        
         

    def homing(self,typ=''):
        """Set up PA
        """
        pass

        #dstatus[13]='-ve limit seen during last move'
        #dstatus[14]='+ve limit seen during last move'
        #while Axes.dstatus[13]!= 1:
        #lire le buffer
        #buf = self.com('R','('+typ+')')
        #buf = buf[1]
        #buf = buf.replace('*','').replace('\r\n','').split('_')
        #print buf
        #if Axes.dstatus[13]== 0:
            #if Axes.dstatus[14]== 1:
                #scom1 = 'D-793600'
                #scom2 = 'G'
                #com   = self.com('D-793600',verbose=True)
                #com   = self.com('GH',verbose=True)
                #com   = self.com('R(PA)',verbose=True)
            #if :
        #else:
            #com   = self.com('R(PA)',verbose=True)

                         
        
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
        self.a  = ['',Axes(1,'x',self.ser,scale=12800),
                      Axes(2,'y',self.ser,scale=22800),
                      Axes(3,'rot',self.ser,scale=2111.1111111111113,typ='r')] #self.a4  = Axes(4,'z',self.ser,typ='r') 





if __name__=="__main__":
    s = Scanner('/dev/ttyUSB0')
    #s = Scanner('/dev/ttyUSB1')
    #sm.fromfile('prog1')
    #sm.fromfile('AY')
    #Sc[1].com('ON')
    #st = sm.com(1,'LIMITS',(0,1,1))
    #st = sm.com(1,'1D-4000')
    #st = sm.com(1,'G')
