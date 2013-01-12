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

    In [1]: H = vc.tran2prop()


Tchannel object is the container for the ray transfer function 

.. ipython::

    In [1]: S.tx.A.info()

    In [1]: plt.figure()

    @savefig rayTF.png width=6in 
    In [2]: H.plot(ix=np.arange(10))


