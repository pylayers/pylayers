import numpy as np
import scipy as sp
import pylab as pl
from scipy import io
import MySQLdb as msql
import csv
from PyUtil import nbint
from matplotlib import pyplot as plt
from PP_Siradel import *
import mplrc.ieee.transaction
from Bsignal import Noise
import MesCEA as mcea
import Layout
from ffnet import ffnet, mlgraph, savenet, loadnet, exportnet
from scikits.learn import cross_val as cross


#~ N=Noise(100,20,-174,0,50).gating(4.589,4.0,'rect')
#~ plt.plot(N.x,N.y)
#~ plt.show()
MATDIR="/media/DONNEES/Pyproject/cir"
TE=1.0# transmitted energy
c=0.3
id_num=0
pos_init=1
pos_end=3
dist_init=4
dist_end=7
vis_init=8
vis_end=11
Etot_init=12
Etot_end=15
Emax_init=16
Emax_end=19
Efirst_init=20
Efirst_end=23
toa_cum_init=24
toa_cum_end=27
toa_max_init=28
toa_max_end=31
toa_th_init=32
toa_th_end=35
toa_stg_init=36
toa_stg_end=39
lqi_init=40
lqi_end=43
taurms_init=44
taurms_end=47
taumoy_init=48
taumoy_end=51
toa_win_init=52
toa_win_end=55
room_init=56

fname='IDCEAMesurenew.txt'
def displot(pt,ph,coltyp):
    """
        plot segments between points
    pt: tail points array (2 x (2*Nseg))
    ph: head points array (2 x (2*Nseg))

    """
    Nseg = shape(pt)[1] 
    pz = empty((2,))
    pn = zeros((2,))
    for i in range(Nseg):
        pz = vstack((pz,pt[:,i],ph[:,i],pn))
    
    m1 = array([0,0,1])
    mask = kron(ones((2,Nseg)),m1)
    pzz = pz[1:,:].T
    vertices = ma.masked_array(pzz,mask)
    plt.plot(vertices[0,:],vertices[1,:],coltyp)


def get_data(fname):
        '''
        read data from txt File
        '''
        
        cs=csv.reader(open(fname),delimiter=' ')
        data=[]
        for row in cs:
            #if int(row[0]) in traj:
                data.append(row)
        #return data
        #data=data[list(traj)]   # identify only some rooms
        # stock data in a list of lists
        for i in arange(0,len(data),1):
                data[i][id_num]=int(data[i][id_num])
                for j in arange(pos_init,dist_end+1):
                        data[i][j]=float(data[i][j])
                for j in arange(vis_init,vis_end+1):
                        data[i][j]=data[i][j]
                for j in arange(Etot_init,toa_win_end+1):
                        data[i][j]=float(data[i][j])
                for j in arange(room_init,room_init+1):
                        data[i][j]=int(data[i][j])
        return data
        
def toa_max(cir,t):
        """
        TOA=time(max(CIR))
        """
        #
        #  ( ) ^2
        #  
        y2 = cir*cir
        maxy2 = max(y2)
        v = pl.find(y2==maxy2)
        toa = t[v[0]]
        return toa
        
def toa_win(cir,t,w):
        """
            Parameters
            ----------
            w : parameter between 0 and 100
            Lei takes w = 9 
        """
        maxbruit = max(cir[0:1000])
        Max = max(cir)
        #print shape(cir)
        nmax = pl.find(cir==Max)
        n   = nmax
        step = Max/1e2
        thre = Max-step
        
        delta =100
        d   = 0
        nint = 0 
        N = np.array([])
        N = np.hstack((N,n))
        # tant delta est plus grande que w% du Max 
        while delta > w*Max/1e2:
            
            u = pl.find(cir > thre)
            hr = pl.find(u>n)
            g = np.delete(u,hr)
        
            if nmax>=6000:
            #set the fenetre=6000*0.005=30ns
                hl = pl.find(g<nmax-6000)
                u = np.delete(g,hl)
            else:
                u = g
            

            n_int = nbint(u) - 1

            if n_int==0:
                thre = thre-step
                d = d+step
            else:
                delta = Max-maxbruit-d-step  
                d = d+step
                n = u[0]
                N = np.hstack((N,n))
                thre = thre-step
            
            if thre<0:
                break   
        if len(N)>=2:     
            nn = N[-2]
        else:
            nn = N[0]
            
        tau = t[nn]
        return tau 

def cdf(x,colsym="",lab="",lw=1,ms=4.0,markevery=1):

        """
        Plot a cumulative density function
        """
        
        x  = sort(x)
        n  = len(x)
        x2 = repeat(x, 2)
        y2 = hstack([0.0, repeat(arange(1,n) / float(n), 2), 1.0])
        plt.plot(x2,y2,colsym,label=lab,linewidth=lw,ms=ms,markevery=markevery)
        #plt.semilogx(x2,y2,colsym,label=lab,linewidth=lw)      
        
def read_mat(Tx, Rx):
    #try:
        mfile=io.loadmat(MATDIR+'/Tx'+str('%03d'%Tx)+'/cir-tx'+str('%03d'%Tx)+'-rx'+str('%03d'%Rx)+'.mat')
        #Tx and Rx positions
        PTx=mfile['Tx']
        PRx=mfile['Rx'+str(Rx)]
        # distance and TOF
        d=np.sqrt(np.sum((PTx - PRx)*(PTx - PRx),axis=0))[0]
        tof=d/0.3
        #Channel response ans timestamp with omni antenna
        ciro=mfile['ciro'+str(Rx)]
        to=mfile['to'+str(Rx)]
        #Channel response ans timestamp with UWB antenna
        cira=mfile['cira'+str(Rx)]
        ta=mfile['ta'+str(Rx)]
        # Add noise to simulated CR
        #~ N=Noise(ta[-1]-ta[0]+(ta[1]-ta[0])/2,1/(ta[1]-ta[0]),-174,NF=-7,R=50).gating(4.589,4.0,'rect')
        #~ cira_old=cira
        if Tx==248 and Rx==352:
            cira=cira.T#+N.y.reshape(shape(cira.T))
        #~ else:
            #~ cira=cira+N.y.reshape(shape(cira))
            
        #~ plt.subplot(311)
        #~ MeasCR=mcea.CEAMesure(Tx-59,Rx-349)
        #~ plt.plot(MeasCR.tdd.ch1.x,MeasCR.tdd.ch1.y,label='meas')
        #~ plt.axis([0,150,-0.00015,0.00015])
        #~ plt.legend()
        #~ plt.subplot(312)
        #~ plt.plot(ta,cira,label='sim')
        #~ plt.axis([0,150,-0.00015,0.00015])
        #~ plt.legend()
        #~ plt.subplot(313)
        #~ plt.plot(ta,cira_old,label='sim without noise')
        #~ plt.axis([0,150,-0.00015,0.00015])
        #~ plt.legend()
        #~ plt.show()
        
        
        
        
        '''plt.plot(ta,10*np.log10(cira**2),'b-')
        plt.xlabel("delay(ns)")
        plt.ylabel("Amplitude (dBm)")
        plt.axis([0,100,-110,-50])
        plt.savefig("cir.jpg",format='jpg',dpi=300)
        plt.close()'''
        # RSSI: Etot
        etot = (ta[1]-ta[0])*sum(cira*pl.conj(cira))
        Et =10*np.log10 (etot)
        #print np.size(Et)
        # RSSI: PL
        PL=10*np.log10(TE-etot)[0]
        # compute TOA
        #toa  = toa_max(cira*cira,ta)[0]-0.7
        toa = toa_win(cira*cira,ta,9)[0]-0.7
            
        return PTx, PRx, d, tof, Et[0], PL, toa
    #~ except:
        #~ print "File does not exist"
        #~ return nan, nan, nan, nan, nan, nan 

#########################################################################################################################
################################################ CREATE THE DATABASE ####################################################
#########################################################################################################################

#~ bb=[]
#~ for o in np.arange(0.0,0.01,0.001):
    #~ aa=[]
    #~ for i in range(1,343):
        #~ print o, i
        #~ aa.append(abs(read_mat(2, i,o)[3]-read_mat(2, i,o)[-1]))
    #~ bb.append(np.mean(aa))
#~ aa=[]
#~ for i in range(1,343):
    #~ print i
    #~ aa.append(read_mat(10, i)[3]-read_mat(10, i)[-1])
    
    
    
HOST='localhost'
USER='root'
PASSWD='sqlsql'
DATABASE='RTsim'

#Connection to DB
db=msql.connect(HOST, USER, PASSWD,DATABASE)
# We'll need a cursor to this database to excute commands
cursor = db.cursor()

#
#Create tables
#

#~ create a new empty table for transmitters (excuted once)
#~ cursor.execute("CREATE TABLE IF NOT EXISTS TxNodes (TxID INT PRIMARY KEY AUTO_INCREMENT, X DOUBLE, Y DOUBLE, Z DOUBLE)")

#~ create a new empty table for Receivers (excuted once)
#~ cursor.execute("CREATE TABLE IF NOT EXISTS RxNodes (RxID INT PRIMARY KEY AUTO_INCREMENT, X DOUBLE, Y DOUBLE, Z DOUBLE)")

#~ create a new empty table for Simulations (excuted once)
#~ cursor.execute("CREATE TABLE IF NOT EXISTS UWBsim (SimID INT PRIMARY KEY AUTO_INCREMENT, TxID INT, RxID INT, RSSI DOUBLE, TOA DOUBLE, TOF DOUBLE, DIST DOUBLE)")

#
# Fill in the tables
#
# Fill in the table TxNodes(executed once)
#~ for Tx in range(1,52):
    #~ print Tx
    #~ M=read_mat(Tx,50)
    #~ X=M[0][0][0]
    #~ Y=M[0][1][0]
    #~ Z=M[0][2][0]
    #~ cursor.execute("INSERT INTO TxNodes(X,Y,Z) VALUES(%s,%s,%s)",(X,Y,Z))

# Fill in the table RxNodes(executed once)
#~ for Rx in range(1,344):
    #~ print Rx
    #~ M=read_mat(2,Rx)
    #~ X=M[1][0][0]
    #~ Y=M[1][1][0]
    #~ Z=M[1][2][0]
    #~ cursor.execute("INSERT INTO RxNodes(X,Y,Z) VALUES(%s,%s,%s)",(X,Y,Z))
    
# Fill in the table UWBsim(executed once)
#~ for Rx in range(1,344):
    #~ for Tx in range(1,52):
        #~ print Rx, Tx
        #~ try:
            #~ M=read_mat(Tx,Rx)
            #~ RSSI=M[4]
            #~ TOA=M[5]
            #~ TOF=M[3]
            #~ dist=M[2]
            #~ cursor.execute("INSERT INTO UWBsim(TxID, RxID, RSSI, TOA, TOF, dist) VALUES(%s,%s,%s,%s,%s,%s)", (Tx, Rx, RSSI, TOA, TOF, dist))
        #~ except:
            #~ cursor.execute("INSERT INTO UWBsim(TxID, RxID, RSSI, TOA, TOF, dist) VALUES(NULL, NULL, NULL, NULL, NULL, NULL)")


#########################################################################################################################
############################################ Simulations Post Processing ################################################
#########################################################################################################################

def RSSI_modeling():
    d=[]
    RSSI=[]
    toa=[]
    for Tx in range(60,362):
        for Rx in range(350,354):
            print 'Tx', Tx, 'Rx', Rx
            M=read_mat(Tx,Rx)
            RSSI.append(M[4])
            d.append(M[2])
            toa.append(M[6])
            '''try:
                M=read_mat(Tx,Rx)
                if M[2]!=0.0 and M[4]>-100.0:
                    RSSI.append(M[4])
                    d.append(M[2])
                    toa.append(M[6])
            except:
                if M[4]<-100:
                    print "simulation file corrupted"
                else:
                    print "No simulation file"'''

    data=get_data(fname)
    
    RSSI_meas=[] # received RSSI
    j=0
    for i in data:
        for k in arange(Etot_init,Etot_end+1):
            RSSI_meas.append(data[j][k])
        j=j+1
               
    dist=[] # distances
    j=0
    for i in data:
        for k in arange(dist_init,dist_end+1):
            dist.append(data[j][k])

        j=j+1
                
    toa_meas=[] # distances
    j=0
    for i in data:
        for k in arange(toa_win_init,toa_win_end+1):
            toa_meas.append(data[j][k])

        j=j+1
                
    PL=np.polyfit(10*np.log10(d),RSSI,1)
    n=-PL[0]
    PL0=PL[1]
    RSSIsh=RSSI+n*10*np.log10(d)-PL0
    shvar=np.var(RSSIsh)
    
    PL_m=np.polyfit(10*np.log10(dist),RSSI_meas,1)
    n_m=-PL_m[0]
    PL0_m=PL_m[1]
    RSSIsh_m=RSSI_meas+n_m*10*np.log10(dist)-PL0_m
    shvar_m=np.var(RSSIsh_m)
    
    PL_meas=[PL0_m, n_m, np.sqrt(shvar_m)]
    PL_sim=[PL0, n, np.sqrt(shvar)]         
    plt.plot(10*np.log10(dist),RSSI_meas, 'b-',label='Measurements',ms=2.0)
    plt.plot(10*np.log10(d),RSSI, 'r*', label="RT Simulations",ms=2.0)
    #plt.plot(RSSI_meas, 'gs',label='Measurements',ms=4.0)
    #plt.plot(RSSI, 'r*', label="RT Simulations",ms=4.0)
    #plt.plot(10*np.log10(sort(dist)),PL_meas[0]-10*PL_meas[1]*np.log10(sort(dist)), 'g-', linewidth=3)
    #plt.plot(10*np.log10(sort(d)),PL_sim[0]-10*PL_sim[1]*np.log10(sort(d)), 'b-', linewidth=3)
    #pt= np.vstack((10*np.log10(asarray(dist)),asarray(RSSI_meas)))
    #ph= np.vstack((10*np.log10(asarray(d)),asarray(RSSI)))
    #displot(pt,ph,'k-')
    #plt.axis([3,14,-80,-40])
    plt.xlabel(r'$10\log_{10}$'+'(d(m))')
    plt.ylabel('RSSI (dB)')
    plt.grid('on')
    plt.legend()
    plt.savefig("PL_sim_meas.pdf", format='pdf', dpi=300)
    plt.close()
    
    cdf(abs(asarray(RSSI_meas)-asarray(RSSI)), "b-", "", 2)
    plt.grid('on')
    plt.xlabel(r'$|RSSI_M/RSSI_{RT}|$'+ ' (dB)')
    plt.ylabel('Probability')
    plt.savefig('CDF_M_vs_RT.pdf', format='pdf', dpi=300)
    plt.close()
    
    cdf(abs(0.3*asarray(toa_meas)-asarray(d)), "b--", "Measurements", 2)
    cdf(abs(0.3*asarray(toa)-asarray(d)), "r-", "RT Simulations", 2)
    
    plt.grid('on')
    plt.xlabel('Ranging error (m)')
    plt.ylabel('Probability')
    plt.axis([0,20,0,1])
    plt.legend()
    plt.savefig('CDF_Ranging_RT.pdf', format='pdf', dpi=300)
    plt.close()
    
    return [PL0,n,np.sqrt(shvar)],[PL0_m,n_m,np.sqrt(shvar_m)], RSSI, d, toa
    

def RSSI_modeling_Rx(Rx):
    d=[]
    RSSI=[]
    toa=[]
    for Tx in range(60,362):
        print 'Tx', Tx, 'Rx', Rx+349
        try:
            M=read_mat(Tx,Rx+349)
            if M[2]!=0.0 and M[4]>-100.0:
                RSSI.append(M[4])
                d.append(M[2])
                toa.append(M[6])
        except:
            if M[4]<-100:
                print "simulation file corrupted"
            else:
                print "No simulation file"

    data=get_data(fname)
    
    RSSI_meas=[] # received RSSI
    j=0
    for i in data:
        RSSI_meas.append(data[j][Etot_init+Rx-1])
        j=j+1
               
    dist=[] # distances
    j=0
    for i in data:
        dist.append(data[j][dist_init+Rx-1])
        j=j+1
        
    toa_meas=[] # distances
    j=0
    for i in data:
        toa_meas.append(data[j][toa_win_init+Rx-1])
        j=j+1   
    
    PL=np.polyfit(10*np.log10(d),RSSI,1)
    n=-PL[0]
    PL0=PL[1]
    RSSIsh=RSSI+n*10*np.log10(d)-PL0
    shvar=np.var(RSSIsh)
    
    PL_m=np.polyfit(10*np.log10(dist),RSSI_meas,1)
    n_m=-PL_m[0]
    PL0_m=PL_m[1]
    RSSIsh_m=RSSI_meas+n_m*10*np.log10(dist)-PL0_m
    shvar_m=np.var(RSSIsh_m)
    
    PL_meas=[PL0_m, n_m, np.sqrt(shvar_m)]
    PL_sim=[PL0, n, np.sqrt(shvar)]         
    #plt.plot(10*np.log10(dist),RSSI_meas, 'gs',label='Measurements',ms=4.0)
    #plt.plot(10*np.log10(d),RSSI, 'r*', label="RT Simulations",ms=4.0)
    plt.plot(np.arange(1,303),RSSI_meas, 'b-*',label='Measurements',ms=3.0,markevery=5,lw=1)
    plt.plot(np.arange(1,303),RSSI, 'r-D', label="RT Simulations",ms=2.0,markevery=5,lw=1)
    #plt.plot(10*np.log10(sort(dist)),PL_meas[0]-10*PL_meas[1]*np.log10(sort(dist)), 'g-', linewidth=3)
    #plt.plot(10*np.log10(sort(d)),PL_sim[0]-10*PL_sim[1]*np.log10(sort(d)), 'b-', linewidth=3)
    #pt= np.vstack((10*np.log10(asarray(dist)),asarray(RSSI_meas)))
    #ph= np.vstack((10*np.log10(asarray(d)),asarray(RSSI)))
    #displot(pt,ph,'k-')
    
    #~ if Rx==1:
        #~ pl.fill([252,259,259,252],[-100,-100,-20,-20],'b',alpha=0.2)
        #~ pl.fill([1,14,14,1],[-100,-100,-20,-20],'b',alpha=0.2)
        #~ pl.fill([101,111,111,101],[-100,-100,-20,-20],'b',alpha=0.2)
        #~ pl.fill([25,33,33,25],[-100,-100,-20,-20],'b',alpha=0.4)
        #~ pl.fill([97,101,101,97],[-100,-100,-20,-20],'b',alpha=0.4)
        #~ pl.fill([188,227,227,188],[-100,-100,-20,-20],'r',alpha=0.2)
    
    
    plt.axis([0,302,-90,-40])
    #plt.xlabel(r'$10\log_{10}$'+'(d(m))')
    plt.annotate('AP'+str(Rx),[5,-45])
    plt.xlabel('Trajectory Index')
    plt.ylabel('RSSI (dB)')
    plt.grid('on')
    if Rx==2:
        plt.legend()
    plt.savefig("PL_sim_meas_Rx"+str(Rx)+".pdf", format='pdf', dpi=300)
    plt.close()
    
    
    plt.plot(np.arange(1,303),abs(0.3*asarray(toa_meas)-asarray(d)), 'b-*',label='Measurements',ms=3.0,markevery=5,lw=1)
    plt.plot(np.arange(1,303),abs(0.3*asarray(toa)-asarray(d)), 'r-D', label="RT Simulations",ms=2.0,markevery=5,lw=1)
    
    #~ if Rx==1:
        #~ pl.fill([252,259,259,252],[-100,-100,-20,-20],'b',alpha=0.2)
        #~ pl.fill([1,14,14,1],[-100,-100,-20,-20],'b',alpha=0.2)
        #~ pl.fill([101,111,111,101],[-100,-100,-20,-20],'b',alpha=0.2)
        #~ pl.fill([25,33,33,25],[-100,-100,-20,-20],'b',alpha=0.4)
        #~ pl.fill([97,101,101,97],[-100,-100,-20,-20],'b',alpha=0.4)
        #~ pl.fill([188,227,227,188],[-100,-100,-20,-20],'r',alpha=0.2)
    
    
    plt.axis([0,302,0,10])
    #plt.xlabel(r'$10\log_{10}$'+'(d(m))')
    plt.annotate('AP'+str(Rx),[5,9])
    plt.xlabel('Trajectory Index')
    plt.ylabel('Ranging error (m)')
    plt.grid('on')
    if Rx==1:
        plt.legend()
    plt.savefig("TOA_sim_meas_Rx"+str(Rx)+".pdf", format='pdf', dpi=300)
    plt.close()
    
    return [PL0,n,np.sqrt(shvar)],[PL0_m,n_m,np.sqrt(shvar_m)], RSSI, d, toa
    
def show_env():
    p0=array([[298626.5107],[2354073.213]])
    pk=array([[298622.3689,298617.4193,298614.2383,298607.736],[2354082.0733,2354088.4029,2354080.9762,2354088.391]])-p0
    plt.rcParams['figure.dpi']=300
    #plt.figure(figsize=(10,8),frameon=False)
    L=Layout.Layout()
    L.loadstr("sircut.str")
    L.loadfur("FurSiradel.ini")
    L.showGs(furniture=True)
    data=get_data(fname)
    X=[]
    Y=[]
    j=0
    for i in data:
        X.append(data[j][pos_init])
        Y.append(data[j][pos_init+1])
        j+=1
    plt.plot(pk[0,:],pk[1,:],'rD',label='AP',markersize=2)
    for k in range(4):
        plt.annotate("AP"+str(k+1),xy=(pk[0,k]+0.1,pk[1,k]+0.1),size=4)
    plt.plot(X,Y,'bo', label='MS',markersize=1,alpha=0.5)
    for i in np.arange(0,len(X),2):
        plt.annotate(str(i+1),xy=(X[i]-0.1,Y[i]+0.1),size=3)
    plt.legend(numpoints=1,loc=3)
    plt.axis('scaled')
    #plt.show()
    plt.savefig('meas_points.pdf',format='pdf', bbox_inches='tight',pad_inches=-0)
    plt.close()
    
def ffnetwork(inputs, layers, outputs):
    conec=mlgraph((inputs, layers, outputs),biases=False)
    net=ffnet(conec)
    '''NX.draw_graphviz(net.graph, prog='dot')
    plt.show()'''
    return net
        
def Learn_position(ldp):
    print "optimizing ..."
    j=4

    data= get_data(fname)
    # data
    data_svm_list=[]
    target_list= []
    for i in arange(0,len(data),1):
        u=asarray(data[i][ldp: ldp+4])
        data_svm_list.append(u)
        v=data[i][pos_init:pos_init+2]
        target_list.append(v)
        
    data_svm=asarray(data_svm_list)
    target=asarray(target_list)
    skf = cross.LeavePOut(len(target), 1)
    #~ skf = cross.KFold(len(target), 10)
    
    net=ffnetwork(j,2,2)
    Perr=[]
    for train, test in skf:
        
        net.train_bfgs(data_svm[train], target[train])
        
        
        for i in arange(0, len(data_svm[test]), 1):
            
            X_new=data_svm[test][i]
            X_pred= net(X_new)
            Perr.append(pl.dist(X_pred,target[test][i].reshape(shape(X_pred)))) 
    print "Neural Networks", mean(Perr)
    cdf(Perr,colsym="b-*",lab="Measurements",lw=2, ms=5.0, markevery=100)
    
    data_svm_list=[]
    target_list= []
    for Tx in range(60,362):
        RSSI=[]
        for Rx in range(350,354):
            M=read_mat(Tx,Rx)
            RSSI.append(M[4])
        u=asarray(RSSI)
        data_svm_list.append(u)
    for i in arange(0,len(data),1):
        v=data[i][pos_init:pos_init+2]
        target_list.append(v)
        
    data_svm=asarray(data_svm_list)
    target=asarray(target_list)
    #~ skf = cross.LeavePOut(len(target), 1)
    #skf = cross.KFold(len(target), 10)
    
    net=ffnetwork(j,2,2)
    Perr=[]
    for train, test in skf:
        
        net.train_bfgs(data_svm[train], target[train])
        
        
        for i in arange(0, len(data_svm[test]), 1):
            
            X_new=data_svm[test][i]
            X_pred= net(X_new)
            Perr.append(pl.dist(X_pred,target[test][i].reshape(shape(X_pred)))) 
    print "Neural Networks", mean(Perr)
    cdf(Perr,colsym="r-D",lab="RT Simulations",lw=2, ms=3.0, markevery=100)
    plt.legend(loc=4)
    plt.grid('on')
    plt.xlabel('Positioning error (m)')
    plt.ylabel('Probability')
    plt.savefig("RT_fingerprint.pdf",format='pdf',dpi=300)
    plt.close()
    
     
