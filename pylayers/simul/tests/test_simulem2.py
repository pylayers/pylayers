#    S = Simul()
#    narg = len(sys.argv)
#    if narg >=2:
#        S.load('where2.ini')
#        itx = int(sys.argv[1])
#        if narg == 4:
#            rx = range(int(sys.argv[2]),int(sys.argv[3]),1)
#        else:
#            srx = range(1,len(S.rx.points.keys())+1,1)
#
#        S.run(itx,srx)
#    else:
#        S.load('where2.ini')
#   try:
#        filename = sys.argv[1]
#        S.filesimul= pyu.getshort(filename)
#        S.load()
        #S.filesimul='where2.simul'
        #S.filesimul='w2-Tx40.simul'
#    S.CIR(arange(Ntx),arange(Nrx),wav,store_level=16+8+4+2+1,alpha=alpha)
#        S.choose()
#        S.run()
#        S.save()
#    S.run()
#    param= {'band': 0.499,
#               'fc': 4.493,
#               'fe': 50,
#               'tresh': 3,
#               'Tw': 30}
#    wav = Waveform('impulse-generic',param)
##    S.save()
#    Ntx = S.tx.N
#    Nrx = len(S.filefield[0])
#    alpha = sqrt(1./30.0)
#    S.CIR(arange(Ntx),arange(Nrx),wav,store_level=16+8+4+2+1,alpha=alpha)
#   #S.indoor.interact()
#   #fitra = S.filetra[0][7]
#   #fitud = S.filetud[0][7]
#   #g0    = GrRay3D()
#   #g0.load(fitra,S.indoor)
#
#   #g1 = GrRay3D.GrRayTud()
#   #g1.load(S.filetud[0][7],S.sl,S.mat)
#   #S =Simulation(filestr='siradel-cut.str2',
#   #        fileslab='siradel.slab',
#   #        filemat='siradel.mat',
#   #        filepalch='def.palch',
#   #        filepatra='def.patra',
#   #        filefreq="def.freq",
#   #        filespaTx="Tx.spa",
#   #        filespaRx="Rx.spa",
#   #        fileantTx="dipole_XZ_NEC2.vsh3",
#   #        fileantRx="dipole_XZ_NEC2.vsh3")
#
#
##  # Propagation channel
##  VC=VectChannel(S,0,0,False)
##
##  # Transmission channel: SC0 (omnis antenna) & SCA (realistic antenna)
##  SCO = VC.vec2scal()
##  SCA = VC.vec2scalA(S.tx.A,S.rx.A)
##
##  # waveform to be transmitted
##  #param= {'band': 3,
##  #       'fc': 5,
##        #       'fe': 100,
##        #       'tresh': 10,
##        #       'Tw': 10}
##        #w = Waveform('impulse-generic',param)
##
##  # Wgamma
##  #f = w.sf.x
##  #ygamma = -1j*0.3/(4*pi*f)
##  #gamm  = bs.FUsignal(f,ygamma)
##  #w.sf=w.st.ftshift()
##  #Wgam  = w.sf*gamm
##  #wgam  = Wgam.ift()
##  #wgam  = Wgam.iftshift()
##
##  # simulated received signal
##  #ro    = SCO.applywavB(Wgam)  # omnis antenna
##  #    = SCA.applywavB(Wgam)  # monocone antenna
##
##  # Siradel measured received signal
##  #filemes = pyu.getlong('SIRADEL_08-07-31_P179.mat','measures')
##  M     = UWBMesure(1,1)
##  r3    = M.tdd.ch3
##
##  # construction of realistic transmitted waveform
##  rr    = M.RAW_DATA.tx
##  tv    = arange(1001)*0.005
##  wbis  = bs.TUsignal(tv,rr)
##  wbis.translate(-0.5)
##  Ewbis  = wbis.energy()
##  wbis.y = wbis.y/sqrt(Ewbis)
##  wbis.y = wbis.y-mean(wbis.y)
##  wbis.Yadd_zeros2l(1900)
##  wbis.Yadd_zeros2r(1100)
##  Wbis     = wbis.ft()
##  f        = Wbis.x
##  ygamma   = -1j*0.3/(4*pi*f)
##  gamm     = bs.FUsignal(f,ygamma)
##  Wgambis  = Wbis*gamm
##
##  robis   = SCO.applywavB(Wgambis)
##
##  rabis    = SCA.applywavB(Wgambis)  # monocone antenna
##  S.CIR(arange(3),arange(82),Wgambis,store_level=16+8+4+2+1,alpha=1.0)
##
##        for k in range(82):
##                I1 = hstack((II,S.CIRa.y))
##  #plot(ro.x,ro.y)
##  #plot(robis.x,robis.y,'r-')
##  #plot(ra.x,ra.y)
##  #plot(rabis.x,rabis.y,'b-')
##  #plot(t,r3,'g-')
##  #show()
##
##
############################
### Create Simulation object
### S = Simulation()
### S.run()
### VCg = VectChannel(S,0)
### alpha = pi*rand()
### beta  = pi*rand()
### gamma = pi*rand()
### Ta = MEulerAngle(alpha,beta,gamma)
### alpha = pi*rand()
### beta  = pi*rand()
### gamma = pi*rand()
### Tb = MEulerAngle(alpha,beta,gamma)
###
### VCl   = VCg2VCl(VCg,Ta,Tb)
###
### SCl   = VCl.vec2scal()
### SCAl  = VCl.vec2scalA(S.tx.A,S.rx.A)
###
### rAl = SCAl.applywavB(Wgam)
###
### figure()
### plot(rAl.x,rAl.y)
### ax=gca()
### ax.ticklabel_format(style='sci',axis='y')
### ax.yaxis.major.formatter.set_powerlimits((0,0))
### xlabel('retard (ns)',fontsize=28)
### ylabel('(V)',fontsize=28)
### xticks(fontsize=24)
### yticks(fontsize=24)
### title('r(t)',fontsize=30)
### show()
###
