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

    In [1]: ip  = EnImpulse([],fc,band,thres,fe)
   
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

.. ipython::

    In [9]: vc = S.VC(1,1)
    
    @savefig doadod.png width=8in
    In [10]: vc.doadod()


*Construction of the ScalChannel*

.. ipython::

    sc = vc.vec2scal()


ScalChannel object contains all the information about the ray transfer function 


S.tx.A.info()


sc.H.plot()


# The antenna can also been taken into account


alpha = 1./sqrt(30)
sca = vc.vec2scalA(S.tx.A,S.rx.A,alpha)
sca.H.plot()


# ## Calculate UWB Channel Impulse Response 


cir = sc.applywavB(S.wav.sfg)


cir.plot()


CIR=cir.esd(mode='unilateral')
CIR.plot()


# This is wrong 
# 
# $\mathbf{Y} = \mathbf{S} \odot \mathbf{W}$

# <codecell>

wgam = S.wav.sfg
Y    = sc.apply(wgam)
tau  = Y.tau0
#print 'tau=',tau
ri   = Y.ft1(500,1)
UH   = Y.symHz(500)
uh   = UH.ifft(1)
UH.plot()
plt.figure()
uh.plot()
#figure()
#ip0.plot()
#figure()
#IP0  = ip0.fft()
#IP0.plot()


ips  = Y.ift(500,1)
t    = ips.x 
ip0  = TUsignal(t,ips.y[0,:])


plot(UH.x,real(UH.y[0,:]),UH.x,imag(UH.y[0,:]))
U0 = FHsignal(UH.x,UH.y[0,:])
u0 = U0.ifft(1)
u1 = ifft(U0.y)
plt.figure()
plot(uh.x,uh.y[0,:]*1000+3)
S.wav.st.plot()


U0.plot()


# # Here is the problem 
# 
# For some reason the Hermitian symmetry forcing is not working here


U1=u0.fft()
g = fft(u1)
plot(abs(g))
plt.figure()
s  = fftshift(u1)
plot(abs(g))
#plt.figure()
#plot(uh.x,uh.y[0,:])
#plot(uh.x,s*50+0.003)

plot(abs(fft(s)),'r')
plot(abs(fft(uh.y[0,:])),'g')


wgam.plot()


S.wav.sf.plot()


print uh.y[0,:]


plot(imag(s))


# Problem $s$ is not real 


u0


plot(real(u0.y))


plot(imag(s))


U0.y


plot(real(U0.y))


U0.y[0]


U0.y[50]


U0.y[-50]


UH.y[0,2]


UH.y[0,-2]


N = len(UH.y)


v1 = UH.y[1:(N-1)/2.]
v2 = UH.y[N:-1:(N-1)/2.]


len(v1)


len(v2)


UH.y[0,-1]


UH.y[0,1]


plot(real(UH.y[0,:]))
plot(imag(UH.y[1,:]))



