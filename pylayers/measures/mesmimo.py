from pylayers.signal.bsignal import *
from pylayers.util.project import *
from pylayers.gis.readvrml import *
import numpy as np 
import matplotlib.pylab as plt 
import matplotlib.animation as animation
import scipy.linalg as la
#
# This class handles the data coming from the MIMO Channel Sounder IETR lab
#
class MIMO(object):
    def __init__(self,
                 _filename='',
                 rep='',
                 Nf=1601,
                 fminGHz=1.8,
                 fmaxGHz=2.2,
                 calibration=True,
                 time=True,
                 Nz=100,
                 Nt=4,
                 Nr=8):
        """
       
        Parameters
        ----------

        filename
        rep 
        fminGHz 
        fmaxGHz
        Nf  
        calibration : Boolean
        Nz : int 
            Number of Zeros
        nT : int 
            (default = 1) 
        
        Notes 
        -----
        
        Data are placed in the directory mesdir + rep directory

        """

        self.freq = np.linspace(fminGHz,fmaxGHz,Nf)

        self.Nt = Nt
        self.Nr = Nr 

        self.Nf = Nf 

        if _filename <> '':
            self.filename = mesdir+rep+_filename
            # load file 
            self.loadraw(self.Nt)
            if calibration:
                self.calibration()
                if time:
                    # reshaping for using ift (todo update ift for MDA !!)
                    Hcal = FUsignal(self.Hcal.x,np.reshape(self.Hcal.y,(Nt*Nr,Nf)))
                    hcal = Hcal.ift(Nz=Nz,ffts=1)
                    shh = hcal.y.shape
                    self.hcal = TUsignal(hcal.x,np.reshape(hcal.y,(Nr,Nt,shh[-1])))

    def __sub__(self,m):
        N = MIMO()
        N.freq = self.freq
        N.Nt = self.Nt
        N.Nr = self.Nr 
        N.Hcal = self.Hcal - m.Hcal
        return(N)

    def loadraw(self,nT=1):
        """ load a MIMO NrxNt raw data sounder file

            The sounder output file is a 2 columns ASCII csv file 
            Module (dB) ;  Angle (Degree)

        """

        fd  = open(self.filename)
        lis = fd.readlines()
        fd.close()
        module  = []
        phasedeg   = []

        for l in lis:
            l.replace('\r\n','')
            g = l.split(';')
            module.append(float(g[0]))
            phasedeg.append(float(g[1]))

        m   = np.array(module)
        phi = np.array(phasedeg)*np.pi/180.
        m   = m.reshape(self.Nr*self.Nt,self.Nf)
        phi = phi.reshape(self.Nr*self.Nt,self.Nf)
        y   = 10**(m/20)*np.exp(1j*phi)
        #
        # Nr x Nt x Nf    (8x4x1601)
        #
        y   = y.reshape(self.Nr,self.Nt,self.Nf)

        self.H = FUsignal(self.freq,y)

    def calibration(self):
        """ Apply calibration files 

        """
        for iR in range(self.Nr):
            for iT in range(self.Nt):
                filename = 'Calib'+str(iT+1)+'x'+str(iR+1)+'.txt'
                C = MIMO(filename,'/calibration/',calibration=False,Nt=self.Nt)
                try:
                    #tc = np.vstack((tc,C.H.y[iR*4+iT,:])) 
                    tc = np.vstack((tc,C.H.y[iR,iT,:])) 
                except:
                    #tc = C.H.y[iR*4+iT,:] 
                    tc = C.H.y[iR,iT,:] 

        tc = tc.reshape(self.Nr,self.Nt,self.Nf)

        self.C = FUsignal(C.freq,tc)

        self.Hcal = self.H/self.C

        del self.H
        del self.C
    
    def calBento2(self,fcGHz=2,duR=0.05,duT=0.05,time=False,taumin=0,taumax=80,Nz=20000):
        """ calculate the Bentosela h function

        """
        # f : m x n x uR x f
        fGHz = self.freq[np.newaxis,np.newaxis,np.newaxis,:]
        # m : m x n x  uR x f
        m = np.arange(self.Nr)[:,np.newaxis,np.newaxis,np.newaxis]
        # uR : m x n x uR x f
        uR = np.arange(-1,1,duR)[np.newaxis,np.newaxis,:,np.newaxis]
        # m : m x n x uR x f
        eR = np.exp(-1j*np.pi*m*uR*fGHz/fcGHz)
        # S : m x n x uR x f
        S  = self.Hcal.y[:,:,np.newaxis,:] * eR  
        #  SR : n x uR x uT x f 
        SR = np.sum(S,axis=0)[:,:,np.newaxis,:]
        # n : n x uR x uT x f
        n = np.arange(self.Nt)[:,np.newaxis,np.newaxis,np.newaxis]
        # uT : n x uR x uT x f
        uT = np.arange(-1,1,duT)[np.newaxis,np.newaxis,:,np.newaxis]
        # eT : n x uR x uT x f
        eT = np.exp(-1j*np.pi*n*uT*fGHz/fcGHz)
        # summation along axix m and n 
        self.h = np.sum(SR*eT,axis=0)

        self.uR = np.arange(-1,1,duR)
        self.uT = np.arange(-1,1,duT)
        NuR = len(self.uR)
        NuT = len(self.uT)
        Nf  = len(self.freq)

        if time: 
            #T = fft.ifft(self.h,axis=2)
            #self.h = abs(fft.fftshift(T,axes=2))
            H = FUsignal(self.freq,np.reshape(self.h,(NuR*NuT,Nf)))
            h = H.ift(Nz=Nz,ffts=1)
            ut = np.where((h.x>taumin) & (h.x<taumax))[0]
            xlim = h.x[ut]
            ylim = h.y[...,ut]
            npts = len(ut)
            self.h = TUsignal(xlim,np.reshape(ylim,(NuR,NuT,npts)))

    def grid(self,M,
             OR=np.array([3.4,0.73]),
             OT=np.array([5.29,6.65]),
             cT=np.array([-0.07,0]),
             cR=np.array([0.07,0])):
        """ Evaluate the data on a grid in the plane 
            
            Parameters 
            ----------

            M : np.array() (Nx x Ny)
            OR : np.array (,2) 
                Origin of receiver [3.4,0.73]
            OT : np.array (,2) 
                Origin of transmitter [5.29,6.65]
            cR : np.array (,2) 
                array receiving vector [0.07,0]
            cT : np.array (,2) 
                array transmitting vector [-0.07,0]

            Notes
            -----

           Updated object members 
            
            self.grid : M  (Nx x Ny x 2) 
            self.gloc : TUsignal (x (,ntau) y (Nx x Ny,ntau) )

        """

        aR = cR[0]/np.sqrt(cR[0]**2+cR[1]**2)
        bR = cR[1]/np.sqrt(cR[0]**2+cR[1]**2)

        aT = cT[0]/np.sqrt(cT[0]**2+cT[1]**2)
        bT = cT[1]/np.sqrt(cT[0]**2+cT[1]**2)
        # mapping 
        uT = (aT*(M[...,0]-OT[0])+bT*(M[...,1]-OT[1]))/np.sqrt((M[...,0]-OT[0])**2+(M[...,1]-OT[1])**2)
        uR = (aR*(M[...,0]-OR[0])+bR*(M[...,1]-OR[1]))/np.sqrt((M[...,0]-OR[0])**2+(M[...,1]-OR[1])**2)
        # sampling in uR and uT 
        uuR = self.uR
        uuT = self.uT
        # index in uR and uT
        iUr=np.array(map(lambda x : np.where(abs(uuR-x)==(abs(uuR-x)).min())[0][0], np.ravel(uR)))
        iUt=np.array(map(lambda x : np.where(abs(uuT-x)==(abs(uuT-x)).min())[0][0], np.ravel(uT)))
        
        self.grid = M
        shM = M.shape
        self.gloc = TUsignal(self.h.x,self.h.y[iUr,iUt,:])
        #shL =  gloc.shape
        #assert(shL[0]==shM[0]*shM[1])
        #self.gloc = np.reshape(gloc,(shM[0],shM[1],shL[1]))



    def capacity(self,trho=np.arange(-40,120,1)):
        """
            log_2(det(I+(Et/N0Nt)HH^{H})
        """
        
        tC = np.zeros((len(trho),len(self.freq)))

        for ir in range(len(trho)):
            rho = 10**(trho[ir]/10)
            for ik in range(len(self.freq)):
                #pdb.set_trace()
                H   = self.Hcal.y[:,:,ik]
                HHH = np.dot(H,np.conj(H).T) 
                IR  = np.eye(self.Nr)
                tC[ir,ik] = np.log(la.det(IR+(rho/self.Nt)*HHH))/np.log(2)
       
        self.C=tC                           


    def mulcplot(self,mode,**kwargs):
        """
        """

        defaults = { 'types' : ['m'],
                   'titles' : np.array([['11','12','13','14'],
                 ['21','22','23','34'],
                 ['31','32','33','34'],
                 ['41','42','43','44'],
                 ['51','52','53','54'],
                 ['61','62','63','64'],
                 ['71','72','73','74'],
                 ['81','82','83','84']]),
                   'ylabels':np.array([['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','','']]),
                   'xlabels':np.array([['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['fGHz','fGHz','fGHz','fGHz']]),
                   'labels':np.array([['calibrated','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','','']])
                 }
       
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        
        if mode=='f':
            fig,ax = self.Hcal.plot(**kwargs)
        else:
            kwargs['xlabels'] = np.array([['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['t(ns)','t(ns)','t(ns)','t(ns)']]),
            fig,ax = self.hcal.plot(**kwargs)

        return(fig,ax)
    
    def showgrid(self,**kwargs):
        """ show the data on a spatial grid 

        Parameters
        ----------
        
        'layout':[],
        's':50,
        'vmin' : 0, 
        'vmax': 0.5,
        'linewidth':0,
        'fig':[],
        'ax':[],
        'save':True,
        'filename':'showgrid1',
        'title':'',
        'save':True,
        'dB':False,
        'OR' : np.array([3.4,0.73]),
        'OT' : np.array([5.29,6.65]),
        'cR' : np.array([0.07,0]),
        'cT' : np.array([-0.07,0]),
        'target' : np.array([]),
        'gating':False,
        'dynamic':30
       

        Notes 
        -----

        This function accvept a Layout as input and allows to display 
        a projection of the spatio-delay volume on a 2D grid. 


        """
        defaults = { 'layout':[],
                    's':50,
                    'vmin' : 0, 
                    'vmax': 0.5,
                    'linewidth':0,
                    'fig':[],
                    'ax':[],
                    'save':True,
                    'filename':'showgrid1',
                    'title':'',
                    'save':True,
                    'dB':False,
                    'OR' : np.array([3.4,0.73]),
                    'OT' : np.array([5.29,6.65]),
                    'cR' : np.array([0.07,0]),
                    'cT' : np.array([-0.07,0]),
                    'target' : np.array([]),
                    'gating':False
                   }
       
       
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        
        OR = kwargs['OR']
        OT = kwargs['OT']
        cR = kwargs['cR']
        cT = kwargs['cT']

        ULAR = OR+np.arange(8)[:,np.newaxis]*cR-3.5*cR
        ULAT = OT+np.arange(4)[:,np.newaxis][::-1]*cT-1.5*cT

        if kwargs['gating']:
            dTM = np.sqrt((self.grid[...,0]-OT[0])**2+(self.grid[...,1]-OT[1])**2)
            dRM = np.sqrt((self.grid[...,0]-OR[0])**2+(self.grid[...,1]-OR[1])**2)
            # dM : Nx,Ny
            dM  = dTM+dRM 
            # dM : ,Nx x Ny
            dM = np.ravel(dM) 
            # 6 sigma = 1/400MHz
            # 6 sigma = 2.5ns 
            # sigma = (2.5/6)
            # alpha = 1/(2 sigma^2) = 2*(2.5)**2/36 = 0.347
            #
            alpha = 0.347
            # Gaussian gate 
            # Laplacian gate 
            # Nx x Ny x Ntau
            self.gate = np.exp(-alpha*(dM[:,np.newaxis]/0.3-self.gloc.x[np.newaxis,:])**2)
            data = self.gloc.y*self.gate
            data = np.sum(abs(data),axis=1)
        else:
            data = np.sum(abs(self.gloc.y),axis=1)

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=(10,10))
            ax  = fig.add_subplot(111)
        else:
            fig=kwargs['fig']
            ax = kwargs['ax']

        if kwargs['dB']:
            data = 20*np.log10(data)
            vmax = data.max()
            # clipping @ vmax - dynamic 
            vmin = vmax-kwargs['dynamic']
        else:
            vmin = data.min()
            vmax = data.max()
        
        scat = ax.scatter(self.grid[...,0],
                               self.grid[...,1],
                               c= data,
                               s=kwargs['s'],
                               vmin=vmin,
                               vmax=vmax,
                               linewidth=kwargs['linewidth'])

        cb = plt.colorbar(scat)
        if kwargs['dB']:
            cb.set_label('Level (dB)')
        else:
            cb.set_label('Linear Level')


        # plot ULAs

        ax.plot(ULAR[:,0],ULAR[:,1],'+b')
        ax.plot(ULAT[:,0],ULAT[:,1],'+g')  
        plt.axis('off')

        # plot target

        if kwargs['target']<>[]:
            target = ax.scatter(kwargs['target'][0],kwargs['target'][1],c='black',s=100)

        # display layout 
        if kwargs['layout'] <> []:
            L = kwargs['layout']
            #fig,ax = L.showG('s',fig=fig,ax=ax,nodes=False)
            L.display['ednodes']=False
            L.display['nodes']=False
            L.display['title']=kwargs['title']
            fig,ax = L.showG('s',fig=fig,ax=ax,nodes=False)

        if kwargs['save']:
            fig.savefig(kwargs['filename']+'.pdf')
            fig.savefig(kwargs['filename']+'.png')

        return fig,ax

    def animgrid(self,**kwargs):
        """
        """

        defaults = { 'layout':[],
                    's':100,
                    'vmin' : 0, 
                    'vmax': 0.5,
                    'linewidth':0,
                    'fig':[],
                    'ax':[],
                    'filename':'animgrid1',
                    'save':True,
                    'abs':True,
                    'title':'',
                   }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        
        if kwargs['fig']==[]:
            fig = plt.figure(figsize=(20,20))
            ax  = fig.add_subplot(111)
        
        if kwargs['layout']<>[]:
            L = kwargs['layout']
            fig,ax = L.showG('s',fig=fig,ax=ax,nodes=False)

        Nframe = self.gloc.y.shape[1]
        if kwargs['abs']:
            scat = ax.scatter(self.grid[...,0],
                               self.grid[...,1],
                               c=abs(self.gloc.y[:,0]),
                               s=kwargs['s'],
                               vmin=kwargs['vmin'],
                               vmax=kwargs['vmax'],
                               linewidth=kwargs['linewidth'])
        else:
            scat = ax.scatter(self.grid[...,0],
                               self.grid[...,1],
                               c=self.gloc.y[:,0],
                               s=kwargs['s'],
                               vmin=kwargs['vmin'],
                               vmax=kwargs['vmax'],
                               linewidth=kwargs['linewidth'])

        title  = ax.text(0.1,0.9,kwargs['title'],transform=ax.transAxes,fontsize=18)
        cb   = plt.colorbar(scat)
        delay_template = '%d : tau = %5.2f (ns) d= %5.2f (m)'
        delay_text  = ax.text(0.1,0.9,'',transform=ax.transAxes,fontsize=18)

        def init():
            delay_text.set_text('')
            if kwargs['abs']:
                scat.set_array(abs(self.gloc.y[:,0]))
            else:
                scat.set_array(self.gloc.y[:,0])
            return scat,delay_text
        
        def animate(i):
            delay_text.set_text(delay_template%(i,self.gloc.x[i],self.gloc.x[i]*0.3))
            if kwargs['abs']:
                scat.set_array(abs(self.gloc.y[:,i]))
            else:
                scat.set_array(abs(self.gloc.y[:,i]))
            return scat,delay_text

        anim = animation.FuncAnimation(fig, 
                                       animate, 
                                       init_func=init, 
                                       frames=Nframe, 
                                       interval=1, 
                                       blit=True)
        if kwargs['save']:
            anim.save(kwargs['filename']+'.mp4', fps=5)
        
        return fig,ax,anim

    def plot(self,frequency=True,phase=False,dB=True,cal=True,fig=[],ax=[],color='k'):
        """

        """
        if fig==[]:
            fig,ax=plt.subplots(8,self.Nt,sharex=True,sharey=True)
        if cal:
            H = self.Hcal
        else:
            H = self.H
        for iR in range(self.Nr):
            for iT in range(self.Nt):
                k = iR*4+iT
                if frequency:
                    if not phase:
                        if dB:
                            #ax[iR,iT].plot(H.x,20*np.log10(abs(H.y[k,:])),color=color) 
                            ax[iR,iT].plot(H.x,20*np.log10(abs(H.y[iR,iT,:])),color=color) 
                        else:
                            #ax[iR,iT].plot(H.x,abs(H.y[k,:]),color='k') 
                            ax[iR,iT].plot(H.x,abs(H.y[iR,iT,:]),color='k') 
                    else:
                        #ax[iR,iT].plot(H.x,np.unwrap(np.angle(H.y[k,:])),color=color) 
                        ax[iR,iT].plot(H.x,np.unwrap(np.angle(H.y[iR,iT,:])),color=color) 
                else:
                        ax[iR,iT].plot(self.h.x,abs(self.h.y[iR,iT,:]),color=color) 
                if (iR==7):         
                    ax[iR,iT].set_xlabel('f (GHz)') 
                ax[iR,iT].set_title(str(iR+1)+'x'+str(iT+1)) 
        return(fig,ax)


#plt.ion()
#TS1 = []
#TS2 = []
#TS3 = []
#TS4 = []
#TS5 = []
#TS1.append(Simo18('S1R0','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R1','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R2','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R3','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R4','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R5','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R6','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R7','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R8','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R9','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R10','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R11','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R12','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R13','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R14','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R15','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R16','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R17','/data/21Novembre2012/S1/'))
#TS1.append(Simo18('S1R18','/data/21Novembre2012/S1/'))
#TS5.append(Simo18('S5R1','/data/21Novembre2012/S5/'))
#TS5.append(Simo18('S5R2','/data/21Novembre2012/S5/'))
#TS5.append(Simo18('S5R3','/data/21Novembre2012/S5/'))
#
#t = TS1[0].h.x
#u = np.nonzero((t>0)&(t<400))[0]
#for k in range(18):
#    try:
#        V = np.dstack((V,TS1[k].h.y))
#    except:
#        V = TS1[k].h.y
#
#plt.figure(figsize=(20,10))
#for k in range(8):
#    plt.subplot(1,8,k+1)
#    plt.imshow(np.log10(abs(V[k,u,:])),origin='lower')
#    plt.axis('auto')
#    plt.axis('off')
#    #plt.title(str(k+1))
#    #plt.xlabel('distance')
#    #plt.ylabel('Time')
#plt.savefig('fig.png')    
##S1R0.Hcal.window('blackman')
##S1R0.Hcal.plotdB(3,phase=False)
##plt.show()
##campdir1 = mesdir+'/data/7Novembre2011'
##campdir2 = mesdir+'/data/14Novembre2011'
##serie1 = me
##sdir+'/data/21Novembre2012/S1'
##m1,p1 = simoload(campdir1+'/Mesure1.txt')
##H1 = calibrate(m1,p1,freq)
###m1,p1 = simoload(serie1+'/S1R0')
###H1 = calibrate(m1,p1,freq)
##m2,p2 = simoload(campdir1+'/Mesure2.txt')
##H2 = calibrate(m2,p2,freq)
##m4,p4 = simoload(campdir1+'/Mesure4.txt')
##H4 = calibrate(m4,p4,freq)
##m5,p5 = simoload(campdir1+'/Mesure5.txt')
##H5 = calibrate(m5,p5,freq)
###H1.plotdB()
###m1,p1 = load('14Novembre2011/Domo1')
###H1 = calibration(m1,p1,freq)
###m2,p2 = load('14Novembre2011/Domo2')
###H2 = calibration(m2,p2,freq)
##
### <codecell>
##
###plt.figure(figsize=(15,10))
###plt.subplot(211)
##H1.plotdB(3,phase=False)
###plt.axis([1.8,2.2,-70,-30])
###plt.title('Door opened')
###plt.subplot(212)
##H2.plotdB(3,phase=False)
###plt.axis([1.8,2.2,-70,-30])
###plt.title('Door closed')
###plt.savefig('fig1.png')
###
#### <codecell>
###
###iant = 5 
###plt.plot(H1.x,abs(H1.y[iant,:]),'b')
###plt.plot(H2.x,abs(H2.y[iant,:]),'r')
####plot(H3.x,abs(H3.y[iant,:]),'g')
####plot(H4.x,abs(H4.y[iant,:]),'k')
####plot(H5.x,abs(H5.y[iant,:]),'c')
###plt.savefig('fig26.png')
###
#### <codecell>
###
####:titre: Fenetrage du signal frequentiel 
##H1.window('blackman')
##H2.window('blackman')
##h2 = H2.ift(Nz=1000,ffts=1)
###
#### <markdowncell>
###
####
#### There is here an interesting phenomena. The caracterization we are doing can almost be considered as a waveguide
#### and on the sensor from the middle, while the door is closed the interference between the different paths is producing
#### an increase of the field.
###
###
##dist = False
##dB   = False
##tmin = 0
##tmax = 60
##plt.figure()
##plt.title('Open/closed door')
##tsp =[]
##for ix in range(8):
##    tsp.append(plt.subplot(2,4,ix+1))
##    h1.plot(iy=ix,tmin=tmin,tmax=tmax,dB=dB,dist=dist,col='b',sharex=tsp[0])
##    h2.plot(iy=ix,tmin=tmin,tmax=tmax,dB=dB,dist=dist,col='r',sharex=tsp[0])
##    plt.axis([tmin,tmax,-0.04,0.04])
##    plt.title('Rx'+str(ix+1))
###plt.savefig('fig10.png')   
###
###
###dist = False
###dB   = False
###tmin = 0
###tmax = 300
###plt.figure()
###plt.title('Open/closed door')
###for ix in range(8):
###    plt.subplot(2,4,ix+1)
###    h1.plot(iy=ix,tmin=tmin,tmax=tmax,dB=dB,dist=False,col='b')
###    h2.plot(iy=ix,tmin=tmin,tmax=tmax,dB=dB,dist=False,col='r')
###    plt.axis([tmin,tmax,-0.04,0.04])
###    plt.title('Rx'+str(ix+1))
###plt.savefig('fig29.png')   
###
###
####h1.gating(14,30)
####h2.gating(14,30)
###
###
###h1.plot(iy=ix,tmin=0,tmax=80,dB=dB,dist=False,col='b')
###h2.plot(iy=ix,tmin=0,tmax=80,dB=dB,dist=False,col='r')    
###
###
###V1=h1.fft()
###V2=h2.fft()
###
###
###ix = 7
###V1.plot(iy=[ix],dB=True,phase=True)
###V2.plot(iy=[ix],dB=True,phase=True)
###
###
###
###
##PDP1 = TUsignal(h1.x,np.sum(np.abs(h1.y)/8,axis=0))
###PDP2 = TUsignal(h2.x,np.sum(np.abs(h2.y)/8,axis=0))
###
###
###PDP1.plot(tmin=0,tmax=450,dB=True,logx=False)
###PDP2.plot(tmin=0,tmax=450,dB=True,col='b',logx=False)
###
###
