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


.. plot::
    :include-source:

    from PyLayers.Signal.Bsignal import *
    from PyLayers.Simul.SimulEM import *
    from matplotlib.pylab import *
    fc     = 4 
    band   = 2
    thresh = 10
    fe     = 100 
    ip     = EnImpulse([],fc,band,thresh,fe)
    ip.plot()
    show()


*Verification of energy normalization in both domains*


    E1  = sum(ip.y*ip.y)*ip.dx()
    print "Integration in time",E1


P = ip.esd()
E2 = sum(P.y)*P.dx()
print "Integration in frequency domain ",E2


# ## Calcul of UWB channel impulse response


S= Simul()
S.load('where2.ini')


st = S.wav.st
sf = S.wav.sf
S.wav.info()


st.plot()
figure()
sf.plot()

# <markdowncell>

*Construction of the VectChannel*



vc = S.VC(1,1)


vc.doadod()


*Construction of the ScalChannel*

# <codecell>

sc = vc.vec2scal()

# <markdowncell>

# The ScalChannel object contain all the information about the ray transfer function 

# <codecell>

S.tx.A.info()

# <codecell>

sc.H.plot()

# <markdowncell>

# The antenna can also been taken into account

# <codecell>

alpha = 1./sqrt(30)
sca = vc.vec2scalA(S.tx.A,S.rx.A,alpha)
sca.H.plot()

# <markdowncell>

# ## Calculate UWB Channel Impulse Response 

# <codecell>

cir = sc.applywavB(S.wav.sfg)

# <codecell>

cir.plot()

# <codecell>

CIR=cir.esd(mode='unilateral')
CIR.plot()

# <markdowncell>

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

# <codecell>

ips  = Y.ift(500,1)
t    = ips.x 
ip0  = TUsignal(t,ips.y[0,:])

# <codecell>

plot(UH.x,real(UH.y[0,:]),UH.x,imag(UH.y[0,:]))
U0 = FHsignal(UH.x,UH.y[0,:])
u0 = U0.ifft(1)
u1 = ifft(U0.y)
plt.figure()
plot(uh.x,uh.y[0,:]*1000+3)
S.wav.st.plot()

# <codecell>

U0.plot()

# <markdowncell>

# # Here is the problem 
# 
# For some reason the Hermitian symmetry forcing is not working here

# <codecell>

U1=u0.fft()
g = fft(u1)
plot(abs(g))
plt.figure()
s  = fftshift(u1)
plot(abs(g))
#plt.figure()
#plot(uh.x,uh.y[0,:])
#plot(uh.x,s*50+0.003)

# <codecell>

plot(abs(fft(s)),'r')
plot(abs(fft(uh.y[0,:])),'g')

# <codecell>

wgam.plot()

# <codecell>

S.wav.sf.plot()

# <codecell>

print uh.y[0,:]

# <codecell>

plot(imag(s))

# <markdowncell>

# Problem $s$ is not real 

# <codecell>

u0

# <codecell>

plot(real(u0.y))

# <codecell>

plot(imag(s))

# <codecell>

U0.y

# <codecell>

plot(real(U0.y))

# <codecell>

U0.y[0]

# <codecell>

U0.y[50]

# <codecell>

U0.y[-50]

# <codecell>

UH.y[0,2]

# <codecell>

UH.y[0,-2]

# <codecell>

N = len(UH.y)

# <codecell>

v1 = UH.y[1:(N-1)/2.]
v2 = UH.y[N:-1:(N-1)/2.]

# <codecell>


len(v1)

# <codecell>

len(v2)

# <codecell>

UH.y[0,-1]

# <codecell>

UH.y[0,1]

# <codecell>

plot(real(UH.y[0,:]))
plot(imag(UH.y[1,:]))

# <codecell>


