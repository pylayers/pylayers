Bsignal
=======

Generation of an Impulse of normalized energy 
---------------------------------------------

 
An energy normalized impulse can be defined as follows
 
.. math::

     p(t)= \frac{\sqrt{2\sqrt{2}}}{\tau\sqrt{\pi}} \cos(2\pi f_c t) e^{-(\frac{t}{\tau})^2}

with

.. math::

      \tau = \frac{2}{B\pi}\sqrt{\frac{\gamma_{dB}\ln{10}}{20}}
 
where :math:`B` is the bandwidth defined at :math:`\gamma_{dB}` and
:math:`f_c` is the central frequency of the pulse.



.. ipython::

    In [1]: from pylayers.signal.bsignal import *

    In [1]: fc = 4 
    
    In [1]: band   = 2

    In [1]: thresh = 10
    
    In [1]: fe     = 100 

    In [1]: ip  = EnImpulse([],fc,band,thresh,fe)
   
    @savefig simple_impulse.png width=5in
    In [1]: ip.plot()

*Verification of energy normalization in both domains*

.. ipython::

    In [1]: Et  = sum(ip.y*ip.y)*ip.dx()

    In [1]: print "Time integration",Et

    In [1]: P = ip.esd()
   
    In [1]: Ef = sum(P.y)*P.dx()

    In [1]: print "Frequency integration",Ef

A simulation file contains the description of an applied waveform. 

.. ipython::

    In [1]: from pylayers.simul.simulem import *

    In [2]: from matplotlib.pylab import *

    In [3]: S = Simul()

    In [4]: S.load('default.ini')

    In [5]: st = S.wav.st

    In [6]: sf = S.wav.sf

    In [7]: S.wav.info()

    @savefig default_wav.png width=5in 
    In [8]: st.plot()


*Construction of the VectChannel*

:ref:`pylayers.antprop.channel.VectChannel`

.. ipython::

    In [9]: vc = S.VC(1,1)
    
    @savefig doadod.png width=8in
    In [10]: vc.doadod()


*Construction of the ScalChannel*

.. ipython::

    In [1]: sc = vc.vec2scal()


ScalChannel object is the container for the ray transfer function 

.. ipython::

    In [1]: S.tx.A.info()

    In [1]: plt.figure()

    @savefig rayTF.png width=6in 
    In [2]: sc.H.plot(ix=np.arange(10))


The antenna is taken into account as follows

.. ipython::

    In [1]: alpha = 1./sqrt(30)

    In [1]: sca = vc.vec2scalA(S.tx.A,S.rx.A,alpha)
    
    In [1]: plt.figure()

    @savefig scalch.png width=6in
    In [1]: sca.H.plot(ix=arange(10))


To evaluate the UWB Channel Impulse Response (CIR), the simulation waveform 
is applied to the ScalChannel



.. ipython::
    
    In [1]: cir = sc.applywavB(S.wav.sfg)
    
    In [1]: plt.figure()

    @savefig cir.png width=6in 
    In [1]: cir.plot()

.. ipython::

    In [1]: CIR = cir.esd(mode='unilateral')
    
    In [1]: plt.figure()
    
    @savefig CIR.png width=6in
    In [1]: CIR.plot(phase=False,dB=True)


    

.. math:: 

    \mathbf{Y} = \mathbf{S} \odot \mathbf{W}


.. ipython::
    
    In [1]: wgam = S.wav.sfg
    
    In [1]: Y    = sc.apply(wgam)

    In [1]: tau  = Y.tau0

    In [1]: print 'tau=',tau

:math:`r_i(t)` is obtained through the :ref:`pylayers.signal.bsignal.FUDsignal.ft1` function 

.. ipython::

    In [1]: ri   = Y.ft1(500,1)

    In [1]: UH   = Y.symHz(500)

    In [1]: uh   = UH.ifft(1)
    
    In [1]: plt.figure()

    @savefig figUH.png width=6in 
    In [1]: UH.plot()

    In [1]: plt.figure()

    @savefig figuh.png width=6in 
    In [1]: uh.plot()
    
    In [1]: ips  = Y.ift(500,1)
    
    In [1]: t    = ips.x 
    
    In [1]: ip0  = TUsignal(t,ips.y[0,:])

    In [1]: plt.figure()

    @savefig figip0.png width=6in 
    In [1]: plot(UH.x,real(UH.y[0,:]),UH.x,imag(UH.y[0,:]))
    
    In [1]: U0 = FHsignal(UH.x,UH.y[0,:])

    In [1]: u0 = U0.ifft(1)

    In [1]: u1 = ifft(U0.y)

    In [1]: plt.figure()

    @savefig figuh0.png width=6in 
    In [1]: plot(uh.x,uh.y[0,:]*1000+3)

    In [1]: S.wav.st.plot()

    In [1]: plt.figure()


    @savefig figU0.png width=6in 
    In [1]: U0.plot()


Here is the problem 
For some reason the Hermitian symmetry forcing is not working here

.. ipython::

    In [1]: U1 = u0.fft()

    In [1]: g = fft(u1)
   
    @savefig figabsg.png width=6in 
    In [1]: plot(abs(g))
    
    In [1]: plt.figure()
    
    In [1]: s  = fftshift(u1)
    
    @savefig figabffts.png width=6in 
    In [1]: plot(abs(fft(s)),'r')
    
    In [1]: plt.figure()

    @savefig figabuhy.png width=6in 
    In [1]: plot(abs(fft(uh.y[0,:])),'g')

    In [1]: wgam.plot()

    In [1]: S.wav.sf.plot()

    In [1]: print uh.y[0,:]

    In [1]: plot(imag(s))


Problem :math:`s` is not real 


.. ipython::

    In [1]: plot(real(u0.y))

    In [1]: plot(imag(s))

    In [1]: U0.y

    In [1]: plot(real(U0.y))

    In [1]: U0.y[0]

    In [1]: U0.y[50]

    In [1]: U0.y[-50]

    In [1]: UH.y[0,2]

    In [1]: UH.y[0,-2]

    In [1]: N = len(UH.y)

    In [1]: v1 = UH.y[1:(N-1)/2.]

    In [1]: v2 = UH.y[N:-1:(N-1)/2.]

