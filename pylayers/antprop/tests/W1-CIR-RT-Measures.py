#!/usr/bin/python
# -*- coding: latin1 -*-
import pdb
import time 
import tkFileDialog
import datetime
from pylayers.simul.simulem import *
from pylayers.measures.mesuwb import *
import matplotlib.pyplot as plt 

def showRT(itx,itxS,Nray=1,fig=[],ax=[]):
	if ax==[]:
        fig = plt.figure()
		ax  = fig.add_subplot('511')
	
	S.L.display['Node']=False
    S.L.display['Scaled']=False
	S.L.display['Thin']=False
	S.show(fig,ax)
	gr0  = GrRay3D()
	gr1  = GrRay3D()
	gr2  = GrRay3D()
	gr3  = GrRay3D()
	gr0.load(S.filetra[0][itx],S.indoor)
	Nray0 = min(Nray,gr0.n)
	gr0.show(ax,arange(0,Nray0,1),col='r')
	gr1.load(S.filetra[1][itx],S.indoor)
	Nray1 = min(Nray,gr1.n)
	gr1.show(ax,arange(0,Nray1,1),col='b')
	gr2.load(S.filetra[2][itx],S.indoor)
	Nray2 = min(Nray,gr2.n)
	gr2.show(ax,arange(0,Nray2,1),col='k')
	gr3.load(S.filetra[3][itx],S.indoor)
	Nray3 = min(Nray,gr3.n)
	gr3.show(ax,arange(0,Nray3,1),col='m')
	chaine = 'Tx '+str(itxS)+ ' Left : Measure - Right : UWB Ray Tracing Simulation '
	title(chaine)
	

def irx2irxs(S,Rx):
	"""
	irx2irxs(S,Tx,irx) : Convert index point of the simulation into Siradel index point 
	
	Warning : There is a permutation between tx and rx for RT and measurement 

		  Tx measurements == Rx Ray Tracing 

	"""	
	npos = S.rx.N
	pt   = S.rx.position 
	nRx  = shape(Rx)[0]
	dico = {}
	for k1 in range(npos):
		#print pt[0:2:,k1]
		for k2 in range(nRx):
			d = pt[0:2,k1]-Rx[k2,0:2]
			u = dot(d,d)
			if u < 1e-5:
				dico[k1]=k2 
			#print Rx[k2,0:2]
        return(dico)		
	

def figRTvsMes(S,dM,irx,irxS,Nray=30,antenna=True,xmin=0,xmax=100):
	"""
	figRTvMes(S,dM,irx,irxS,Nray=30,antenna=True,xmin=0,xmax=100)

	S    : RT simulation 
	dM   : Measure dictionnary 
	irx  : Transmitter index (start on 0 )
	irxS : Transmitter index (Siradel convention start on 1 )  

	"""
	fig=plt.figure()
	ax=fig.add_subplot(511)
	showRT(irx,irxS,Nray,fig,ax)
	TRT     = Tdd()
	if antenna:
		S.CIRa[0][irx].zlr(0,200)
		S.CIRa[1][irx].zlr(0,200)
		S.CIRa[2][irx].zlr(0,200)
		S.CIRa[3][irx].zlr(0,200)
		TRT.ch1 = S.CIRa[0][irx]
		TRT.ch2 = S.CIRa[1][irx]
		TRT.ch3 = S.CIRa[2][irx]
		TRT.ch4 = S.CIRa[3][irx]
	else:
		S.CIRo[0][irx].zlr(0,200)
		S.CIRo[1][irx].zlr(0,200)
		S.CIRo[2][irx].zlr(0,200)
		S.CIRo[3][irx].zlr(0,200)
		TRT.ch1 = S.CIRo[0][irx]
		TRT.ch2 = S.CIRo[1][irx]
		TRT.ch3 = S.CIRo[2][irx]
		TRT.ch4 = S.CIRo[3][irx]

	dM[irxS].show(display=False,col=['r','b','k','m'],C=1,NC=2,xmin=xmin,xmax=xmax)

	TRT.show(display=False,col=['r','b','k','m'],C=2,NC=2,xmin=xmin,xmax=xmax)
	
	sz   = fig.get_size_inches()
	fig.set_size_inches(sz*2.3)
	ext1 = '.pdf'
	ext2 = '.tex'
	ext3 = '.png'
       	filename='Txa'+str(irx)
	fig.savefig('./where/'+filename+ext3,orientation='portrait')
	fig.clf()
	plt.close(fig)
	#plt.clf()


if (__name__=="__main__"):
	#
	# Select trajectory 
	#
	traject_BD  = array([0,1,2,3,6,7,8,9,10,11,12,19,20,41,44,47,50,53,56,59,62,65,68,71,74,77,112,115,118,120,122,123,124,125,126,127,128,129,130,134,139,142,147,149,161,164,169,172,177,180,185,187,190,191,194,195,198,199,202,204,205,208,209,227,228,231,233,235,237,239,241,243,262,263,264,265,266,267,268,269,270,285])

	traject_IDs = array([1, 2, 3, 5, 8, 9, 10, 11, 12, 13, 14, 21, 22, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 123, 126, 129, 134, 137, 138, 139, 140, 141, 142, 143, 144, 145, 163, 168, 171, 176, 179, 201, 204, 209, 212, 217, 220, 228, 230, 233, 234, 237, 238, 241, 242, 245, 247, 248, 251, 252, 297, 298, 301, 303, 305, 307, 309, 311, 313, 333, 334, 335, 336, 337, 338, 339, 340, 341, 356])
	
    select = arange(len(traject_IDs))
	Nptraj = len(traject_IDs[select])

	zh1    = 1.2*ones(Nptraj).reshape(Nptraj,1)
	zh2    = 1.5*ones(Nptraj).reshape(Nptraj,1)

	pt_h1  = hstack((Rx[traject_IDs[select],:],zh1))
	#pt_h2  = hstack((Rx[traject_IDs[select],:],zh2))
	#pt_h   = vstack((pt_h1,	pt_h2))

	#
	# load simulation with siradel files
	#S.load('simul4')
	Tx,Rx  = ptSiradel()
	# 
	# Read measurement 
	#
	dM = {}
	for k in traject_IDs[select]:
		M=CEAmesure(k,1)
		dM[k]=M
	wCEA = TUsignal(dM[k].tdd.tx.x[0:1000],dM[k].tdd.tx.y[0:1000])
	wCEA.zlr(-5,5)
	WCEA=wCEA.ftshift()
	#fichsimulin = 'simul4'
	#root = Tk()
	#root.withdraw()
	FD = tkFileDialog
	fichsimulin = FD.askopenfilename(filetypes = [("Fichiers simul ","*.simul"),
			("All", "*")],
			title="Please choose a simulation file",
			initialdir=simuldir # parent = root
			)
	_fichsimulin = getshort(fichsimulin)
	S = Simulation(_fichsimulin)
	if S.progress<4:
		isauv_new=True
	else:
		isauv_new=False

	#S.run(4)
	if isauv_new:
		stime  = datetime.datetime.now().strftime("--%d-%m-%Y-%Hh%M")
		chaine1 = _fichsimulin.split('--')
		chaine2 = chaine1[0].split('.')
		racine  = chaine2[0]
		fichsimulout  = racine+stime+'.simul'
		S.save(fichsimulout)

	## waveform to be transmitted
	param= {'band': 3,
	      'fc': 5,
	      'fe': 100,
	       'tresh': 10,
	       'Tw': 10}
	w = Waveform('impulse-generic',param)
	#
	# Wgamma		
	#
	f = w.sf.x
	ygamma = -1j*0.3/(4*pi*f)
	gamm  = FUsignal(f,ygamma)
	w.sf  = w.st.ftshift()
	#Wgam  = w.sf*gamm
	Wgam  = WCEA*gamm
	#showRT(0,10)
	mask=1+2+4+8+16
	irx   = range(S.rx.N)
	alpha = sqrt(1.0/30.)
	itx   = arange(4)
	S.CIR(itx,irx,Wgam,store_level=mask,alpha=alpha)
	#
	# Mise au format Tdd ( Mesure SIRADEL) 
	#
	TTRTo= []
	TTRTa= []
	dirx = irx2irxs(S,Tx)
	for k in irx:
		if (mask&8)==8:
			TRTo     = Tdd()
			if (S.CIRo[0][k]!=[]):
				S.CIRo[0][k].zlr(0,200)
				S.CIRo[1][k].zlr(0,200)
				S.CIRo[2][k].zlr(0,200)
				S.CIRo[3][k].zlr(0,200)
				TRTo.ch1 = S.CIRo[0][k]
				TRTo.ch2 = S.CIRo[1][k]
				TRTo.ch3 = S.CIRo[2][k]
				TRTo.ch4 = S.CIRo[3][k]
				TTRTo.append(TRTo)
			else:
				TTRTo.append([])
		if (mask&16)==16:
			TRTa     = Tdd()
			if (S.CIRa[0][k]!=[]):
				S.CIRa[0][k].zlr(0,200)
				S.CIRa[1][k].zlr(0,200)
				S.CIRa[2][k].zlr(0,200)
				S.CIRa[3][k].zlr(0,200)
				TRTa.ch1 = S.CIRa[0][k]
				TRTa.ch2 = S.CIRa[1][k]
				TRTa.omp        
				TRTa.ch4 = S.CIRa[3][k]
				TTRTa.append(TRTa)
			else:
				TTRTa.append([])
        for k in dM.keys():
                try:
                        M1ch1=vstack((M1ch1,dM[k].tdd.ch1.y))
                        M1ch2=vstack((M1ch2,dM[k].tdd.ch2.y))
                        M1ch3=vstack((M1ch3,dM[k].tdd.ch3.y))
                        M1ch4=vstack((M1ch4,dM[k].tdd.ch4.y))
                except:
                        M1ch1=dM[k].tdd.ch1.y
                        M1ch2=dM[k].tdd.ch2.y
                        M1ch3=dM[k].tdd.ch3.y
                        M1ch4=dM[k].tdd.ch4.y


        for k in irx:
                S.CIRa[0][k].zpl(0,200)
                S.CIRa[1][k].zpl(0,200)
                S.CIRa[2][k].zpl(0,200)
                S.CIRa[3][k].zpl(0,200)
                S.CIRa[0][k].zpl(0,200)
                S.CIRa[1][k].zpl(0,200)
                S.CIRa[2][k].zpl(0,200)
                S.CIRa[3][k].zpl(0,200)


	#figRTvsMes(S,dM,irx,dirx[irx],Nray=30,antenna=True,xmin=0,xmax=100)
	#for irx in range(S.rx.N): 
	#	print irx,dirx[irx]
	#	figRTvsMes(S,dM,irx,dirx[irx],antenna=True)
	#fig=plt.figure()
	#ax=fig.add_subplot(511)
	#showRT(1,30,fig,ax)
	#"dM[179].show(display=False,col=['r','b','k','m'],C=1,NC=2,xmin=0,xmax=100)
	#TTRTo[1].show(display=False,col=['r','b','k','m'],C=2,NC=2,xmin=0,xmax=100)
	#r0.show(ax)
	##
	##  Propagation channel
	#VC = VectChannel(S,0,0,transpose)
	##
	#
	## Transmission channel: SC0 (omnis antenna) & SCA (realistic antenna)
	#
	#
	#wgam  = Wgam.ift()
	#wgam  = Wgam.iftshift()
	#
	## simulated received signal
	#
	# Siradel measured received signal 
	#filemes = getlong('SIRADEL_08-07-31_P179.mat','measures')
	#M       = CEAMesure(filemes)
	#t       = 1e9*M.tdd.time
	##r3     = M.tdd.ch3
	#
	## construction of realistic transmitted waveform 	
	##rr    = M.RAW_DATA.tx
	##tv    = arange(1001)*0.005
	##wbis  = TUsignal(tv,rr)
	##wbis.translate(-0.5)
	##Ewbis  = wbis.energy()
	##wbis.y = wbis.y/sqrt(Ewbis)
	##wbis.y = wbis.y-mean(wbis.y)
	##wbis.Yadd_zeros2l(1900)	
	##wbis.Yadd_zeros2r(1100)	
	##Wbis = wbis.ft()
	##f = Wbis.x
	##ygamma = -1j*0.3/(4*pi*f)
	##gamm  = FUsignal(f,ygamma)
	##Wgambis  = Wbis*gamm
	#
	##robis   = SCO.applywavB(Wgambis)
	#
	##rabis    = SCA.applywavB(Wgambis)  # monocone antenna
	##plot(ro.x,ro.y)
	##plot(robis.x,robis.y,'r-')
	##plot(ra.x,ra.y)
	##plot(rabis.x,rabis.y,'b-')
	##plot(t,r3,'g-')
	##show()
	#
	#
	###########################	
	## Create Simulation object 
	##	S = Simulation()
	##	S.run()
	##	VCg = VectChannel(S,0)
	##	alpha = pi*rand()
	##	beta  = pi*rand() 
	##	gamma = pi*rand() 
	##	Ta = MEulerAngle(alpha,beta,gamma)
	##	alpha = pi*rand() 
	##	beta  = pi*rand() 
	##	gamma = pi*rand() 
	##	Tb = MEulerAngle(alpha,beta,gamma)
	##
	##	VCl   = VCg2VCl(VCg,Ta,Tb)
	##
	##	SCl   = VCl.vec2scal()
	##	SCAl  = VCl.vec2scalA(S.tx.A,S.rx.A)
	##	
	##	rAl = SCAl.applywavB(Wgam)
	##
	##	figure()
	##	plot(rAl.x,rAl.y)
	##	ax=gca()
	##	ax.ticklabel_format(style='sci',axis='y')
	##	ax.yaxis.major.formatter.set_powerlimits((0,0))	
	##	xlabel('retard (ns)',fontsize=28)
	##	ylabel('(V)',fontsize=28)
	##	xticks(fontsize=24)
	##	yticks(fontsize=24)
	##	title('r(t)',fontsize=30)
	##	show()
	##
