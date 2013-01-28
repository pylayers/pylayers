import time

FD = tkFileDialog
      filename = FD.askopenfilename(filetypes = [("Fichiers vsh3","*.vsh3"),
      ("All", "*")],
       title="Please choose a vsh3 antenna file",
       initialdir=antdir)
A=Antenna('vsh3',_filename)

FD = tkFileDialog
filename = FD.askopenfilename(filetypes = [("Fichiers trx","*.trx"),
     ("All", "*")],
title="Please choose a trx antenna file",
    initialdir=antdir)
_filename=os.path.split(filename)[1]
#_filetrx='dipole_XZ_NEC2.trx'
#A=Antenna(_filetrx,'bcp')


toc1 = time.time()
A    = Antenna('trx',_filename)
tic1 = time.time()
print "time needed for reading  ",A._filename,"  file = ",tic1-toc1
toc2=time.time()

if A.tau!=0:
    print 'Electrical delay needed'
    A.elec_delay()

A.vsh()

print "A.C.s1tos2(-1)"
#A.C.s1tos2(-1)
A.C.s1tos2(30)
tic2=time.time()
print "time needed for computing s2 coeff = ",tic2-toc2

print "A.savevsh2()"
A.savevsh2()

print "ifreq=40" 
ifreq = 40  # A.fa[40]=4GHz
th    = np.kron(A.theta,np.ones(A.Np))
ph    = np.kron(np.ones(A.Nt),A.phi)

Fth2,Fph2 = A.Fsynth2(th,ph)
# calcul d'erreur relative de reconstruction 
Err=A.mse(Fth2,Fph2)
print "Err = ",Err

FTh2 = Fth2.reshape(A.Nf,A.Nt,A.Np)
FPh2 = Fph2.reshape(A.Nf,A.Nt,A.Np)

print "compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'modulus')"
compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'modulus')
print "compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'real')"
compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'real')
print "compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'imag')"
compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'imag')
print "compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'phase')"
compdiag(ifreq,A,A.theta,A.phi,FTh2,FPh2,'phase')
print "A.C.show(type='s1',k=ifreq)"
A.C.show(type='s1',k=ifreq)
print "A.C.show(type='s2',k=50)"
A.C.show(type='s2',k=20)
print "A.show3(k=ifreq,type='Gain')"
A.show3(k=ifreq,type='Gain')

print "A.C.s2tos3(1e-10)"
A.C.s2tos3(1e-20)

A.C.show(type=3)
Fth3,Fph3 = A.Fsynth3(th,ph)
FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)
compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,typ='modulus')
compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,typ='phase')
show()

print "minimisation of s3"

MSE   = array([])
MSEth = array([])
MSEph = array([])
TI    = array([])
K     = array([])

tic       = time.time()
Fth3,Fph3 = A.Fsynth3(th,ph)
toc       = time.time()
Err       = A.mse(Fth3,Fph3,0)
nk        = shape(A.C.Br.s3)[1]

      MSE    = np.hstack((MSE,Err[0]))
       MSEth  = np.hstack((MSEth,Err[1]))
       MSEph  = np.hstack((MSEph,Err[2]))
      TI     = np.hstack((TI,toc-tic))
     K      = np.hstack((K,nk))

      Enc  = A.C.ens3()
     n    = len(Enc)

      pos=0
      #while (pos<n)&(Err[0]<0.1):
lim=0.1
#lim=0.05
      while (pos<n)&(Err[0]<lim):

         Emin=Enc[pos]
             d=A.C.drag3(Emin)

              tic  = time.time()
               Fth3,Fph3 = A.Fsynth3(th,ph)
              toc  = time.time()
               Err  = A.mse(Fth3,Fph3,0)
     
    if Err[0]>=lim:
        i=d[0][0]
                    i3=d[1][0]
                      A.C.put3(i,i3)
                      Fth3,Fph3 = A.Fsynth3(th,ph)
                      Err  = A.mse(Fth3,Fph3,0)

    nk = shape(A.C.Br.s3)[1]

               MSE    = np.hstack((MSE,Err[0]))
               MSEth  = np.hstack((MSEth,Err[1]))
              MSEph  = np.hstack((MSEph,Err[2]))
               TI   = np.hstack((TI,toc-tic))
              K    = np.hstack((K,nk))
               pos=pos+1


       FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
       FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)

      compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,'modulus','english')
       compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,'phase','english')

       figure()
       plot(K,TI,linewidth=3)
       xlabel('index',fontsize=30)
      ylabel('Time (seconds)',fontsize=30)
      #xlabel('index',fontsize=30)
       #ylabel('Time (s)',fontsize=30)
       xticks(fontsize=26)
       yticks(fontsize=26)
       grid(linewidth=2)

       figure()
      plot(K,100*MSEth,'b-',K,100*MSEph,'r-',linewidth=3)
       xlabel('index ',fontsize=30)
#          ylabel('Error (percent)',fontsize=30)
       #xlabel('index',fontsize=30)
       #ylabel('Error (absolute value)',fontsize=30)
       xticks(fontsize=26)
       yticks(fontsize=26)
       grid(linewidth=2)
      legend((r'${F_{\theta}}$',r'${F_{\phi}}$'),prop=FontProperties(size=30))

show()


#### 2eme etape: Creation de la liste A, contenant toutes les antennes intialisees avec les fichiers vsh2 existents
A=[]
dmin =20  # longeur min du plan de masse simule (mm)
dmax =112
pas  =4
D=np.arange(dmin,dmax,pas)
for index in D:
    _filename='monopole_GroundPlane'+str(index)+'mm.vsh2'
       ant=Antenna('vsh2',_filename)
         A.append(ant)


#############################
Etude Complet!!!    
#############################
Attention: probleme au niveau de la minimisation de shape3

toc1=time.time()
A=Antenna('trx',_filename)
tic1=time.time()
print "time needed for reading  ",A._filename,"  file = ",tic1-toc1

toc2=time.time()
if A.tau!=0:
    print 'Electrical delay needed'
    A.elec_delay()
A.vsh()

print "A.C.s1tos2(30)"
A.C.s1tos2(30)
print "A.C.s2tos3(1e-15)"
A.C.s2tos3()
print "minimise s3 and create vhs3 file"
A.minsh3(emax=0.05)
tic2=time.time()

print "time needed for computing s3 coeff = ",tic2-toc2

A.savevsh3()

print "ifreq=0" 
ifreq=0  # A.fa[0]=1GHz
th    = np.kron(A.theta,np.ones(A.Np))
ph    = np.kron(np.ones(A.Nt),A.phi)

Fth3,Fph3 = A.Fsynth3(th,ph)

       FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
       FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)

print "compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,'modulus')"
    compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,'modulus')
print "compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,'phase')"
      compdiag(ifreq,A,A.theta,A.phi,FTh3,FPh3,'phase')

