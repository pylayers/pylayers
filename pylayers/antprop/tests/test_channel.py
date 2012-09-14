###########################
# test Channel, Simulation, Antenna
###########################
from pylayers.antprop.channel import *
import numpy as np 
#
# choix des antennes
#
fileantt = FD.askopenfilename(filetypes = [("Fichiers vsh3","*.vsh3"),
     ("All", "*")],
      title="Please choose an antenna file",
      initialdir=antdir)

fileantr = FD.askopenfilename(filetypes = [("Fichiers vsh3","*.vsh3"),
("All", "*")],
title="Please choose an antenna file",
initialdir=antdir)

#
# create a Simulation object 
#
S = Simulation(fileantTx=_fileantt,fileantRx=_fileantr)

#
# choix du fichier field
#
filefield = FD.askopenfilename(filetypes = [("Fichiers field","*.field"),
("All", "*")],
title="Please choose a field file",
initialdir=tuddir)

filetauk = filefield.replace('.field','.tauk')
filetang = filefield.replace('.field','.tang')
filerang = filefield.replace('.field','.rang')
filefreq = filefield.replace('.field','.freq')
# il faudra recuperer aussi filetra

S.filefield.append(getshort(filefield))
S.filetauk.append(getshort(filetauk))
S.filetang.append(getshort(filetang))
S.filerang.append(getshort(filerang))
#S.filetra.append(getshort(filetra))

fichsimulin = FD.askopenfilename(filetypes = [("Fichiers simul ","*.simul"),
  ("All", "*")],
  title="Please choose a simulation file",
  initialdir=simuldir)
#          parent=root)
_fichsimulin = getshort(fichsimulin)
S = Simulation(_fichsimulin)

#
# Frequency range for Siradel measurements
#
f     = arange(2,11.05,0.05)

#################################
# definition de la forme d'onde a l'emission
#################################
Tw = 10
x  = arange(-Tw,Tw,0.005)
w  = EnImpulse(x,5,3,10)

W = w.ft()
f = W.x
ygamma = -1j*0.3/(4*np.pi*f)
gamm  = FUsignal(f,ygamma)
Wgam  = W*gamm
wgam  = Wgam.ift()

##       #
##       # plot waveform
##       #
##       figure()
##       plot(w.x,w.y,linewidth=3)
##       ylabel('(V)',fontsize=28)
##       xlabel('retard(ns)',fontsize=28)
##       xticks(fontsize=26)
##       yticks(fontsize=26)
##       title('w(t)',fontsize=30)
##       grid(linewidth=2)
##       show()
#
#      ################################
#      # recuperation du canal de propagation simule
#      ################################
#
#      VC   = VectChannel(S,0,False)
#      #
#      # If file from Siradel
#      #
#      #VC   = VectChannel(S,0,True)
#
#      ################################
#      # calcul du canal de transmission
#      ################################
#
#      #
#      # avec antennes omnis
#      #
#      SC      = VC.vec2scal()
#
#      #
#      # avec antennes simulees
#      #
#      VCA     = VC
#      SCA     = VCA.vec2scalA(S.tx.A,S.rx.A)
#
#      ################################
#      # calcul du signal recu
#      ################################
#
#      #
#      # avec antennes omnis
#      #
#      r  = SC.applywavB(Wgam)
#
#      #
#      # avec antennes simulees
#      #
#      rA = SCA.applywavB(Wgam)
#
#      ################################
#      # visualisation de resultats
#      ################################
#
#      figure()
#      # pour monocone on a un gain max de 8.59**2
#      #plot(r.x,r.y*70,'b-')
#      # pour dipole on a un gain max de 1.68**2
#      #plot(r.x,r.y*2.8,'b-')
#
#      #
#      # calcul de la valeur maximale du gain d'antenne
#      #
#
#      gain_At=np.sqrt((abs(SCA.Ftt.y))**2+(abs(SCA.Ftp.y))**2)
#      gain_Ar=np.sqrt((abs(SCA.Frt.y))**2+(abs(SCA.Frp.y))**2)
#
#      G=gain_At.max()*gain_Ar.max()
#
#      figure()
#      plot(r.x,r.y*G,'b-')
#      plot(rA.x,rA.y,'r-')
#      ax=gca()
#      ax.ticklabel_format(style='sci',axis='y')
#      ax.yaxis.major.formatter.set_powerlimits((0,0))
#      xlabel('retard (ns)',fontsize=28)
#      ylabel('(V)',fontsize=28)
#      xticks(fontsize=24)
#      yticks(fontsize=24)
#      title('r(t)',fontsize=30)
#      show()
#
##       #legend(('antennes omnis','antennes simulees'))
##       show()
#
###############################
## test changement d'antennes
###############################
#      #S.tx.gantenna()
#      #S.rx.gantenna()
#
#      #SCA   = VCA.vec2scalA(S.tx.A,S.rx.A)
#      #rA = SCA.applywavB(Wgam)
#
#      #gain_At=np.sqrt((abs(SCA.Ftt.y))**2+(abs(SCA.Ftp.y))**2)
#      #gain_Ar=np.sqrt((abs(SCA.Frt.y))**2+(abs(SCA.Frp.y))**2)
#
#      #G=gain_At.max()*gain_Ar.max()
#
#      #figure()
#      #plot(r.x,r.y*G,'b-')
#      #plot(rA.x,rA.y,'r-')
#      #xlabel('retard (ns)',fontsize=28)
#      #ylabel('(V)',fontsize=28)
#      #xticks(fontsize=24)
#      #yticks(fontsize=24)
#      #title('r(t)',fontsize=30)
#
#      #legend(('ideale','dipole'))
#      #show()
#
#

##############################
# test passage Rg -> Rl
##############################
itx=0
irx=0
VCg = VectChannel(S,itx,irx,False)
alpha = np.pi*rand()
beta  = np.pi*rand()
gamma = np.pi*rand()
Ta = MEulerAngle(alpha,beta,gamma)
alpha = np.pi*rand()
beta  = np.pi*rand()
gamma = np.pi*rand()
Tb = MEulerAngle(alpha,beta,gamma)

VCl=VCg2VCl(VCg,Ta,Tb)

SCl  = VCl.vec2scal()

VClA  = VCl
SClA  = VClA.vec2scalA(S.tx.A,S.rx.A)

rl = SCl.applywavB(Wgam)
rlA = SClA.applywavB(Wgam)

gain_At=np.sqrt((abs(SClA.Ftt.y))**2+(abs(SClA.Ftp.y))**2)
gain_Ar=np.sqrt((abs(SClA.Frt.y))**2+(abs(SClA.Frp.y))**2)
G=gain_At.max()*gain_Ar.max()

figure()
plot(rl.x,rl.y*G,'b-')
plot(rlA.x,rlA.y,'r-')
ax=gca()
ax.ticklabel_format(style='sci',axis='y')
ax.yaxis.major.formatter.set_powerlimits((0,0))
xlabel('retard (ns)',fontsize=28)
ylabel('(V)',fontsize=28)
xticks(fontsize=24)
yticks(fontsize=24)
title('r(t)',fontsize=30)
show()

