```python
>>> %pylab inline
```

## Synthesis of Ultra Wide Band Waveforms

Once the propagation channel has been evaluated. This is done in the `pylayers.antprop.channel` module.
The received signal is evaluated in applying a convolution product of each ray tranfer function with a specific
IR-UWB waveform. The necessary modules are

+ `pylayers.signal.bsignal`.
+ `pylayers.signal.waveform`
+ `pylayers.signal.channel`

The module `pylayers.simul.simulem` is for definition of electromagnetic simulation.

```python
>>> from pylayers.signal.bsignal import *
>>> from pylayers.signal.waveform import *
>>> from pylayers.antprop.channel import *
>>> from pylayers.simul.simulem import *
```

### Generation of an Impulse of normalized energy

One possible manner to define an energy normalized short UWB impulse is as follows using  `bsignal.EnImpulse` function.

The default waveform is a gaussian windowing of a sine wave of frequency $f_c$. The normalization term depends on the exponential scaling factor $\tau$.

$p(t)= \frac{\sqrt{2\sqrt{2}}}{\tau\sqrt{\pi}} \cos(2\pi f_c t) e^{-(\frac{t}{\tau})^2}$

$\tau = \frac{2}{B\pi}\sqrt{\frac{\gamma_{dB}\ln{10}}{20}}$

where $B$ is the desired bandwidth defined at $\gamma_{dB}$ below the spectrum maximum and $f_c$ is the central frequency of the pulse.

```python
>>> fc     = 4
>>> band   = 2
>>> thresh = 10
>>> fe     = 100
>>> ip     = EnImpulse([],fc,band,thresh,fe)
```

```python
>>> ip.info()
```

## Verification of energy normalization in both domains

```python
>>> E1= sum(ip.y*ip.y)*ip.dx()
>>> print "Integration in time",E1
```

```python
>>> P = ip.esd()
>>> E2 = sum(P.y)*P.dx()
>>> print "Integration in frequency domain ",E2
```

## Calculation of UWB channel impulse response

We choose to load a simple floor plan.

```python
>>> S = Simul()
>>> S.L = Layout('defstr3.ini')
```

A simulation object has an `info` method providing a summary of simulation informations.

```python
>>> st = S.wav.st
>>> sf = S.wav.sf
>>> S.wav.info()
```

The waveform associated with the simulation object is

```python
>>> S.wav
```

```python
>>> S.wav.show()
```

Above the waveform is a generic UWB waveform. The interested user can add easyly any other mathematical expression of UWB waveform for investigation on pulse waveform modulation for example. The waveform can also comes from measurement. For now there are two version of this waveform which has been used during the M1 measurement campaign. One is not compensated `W1compensate` for an extra short delay which can introduse a bias when interpretating the observed delay in terms of distance. The non compensated version is `W1offset` from the time origin about 0.7 ns.

The waveform class should grow for incorporating more waveforms, especially waveforms compliants with the current IEEE 802.15.4a and IEEE 802.15.6 standards.

```python
>>> wavmeasured = Waveform(typ='W1compensate')
>>> wavmeasured.show()
```

```python
>>> wavmeasured = Waveform(typ='W1offset')
>>> wavmeasured.show()
```

Here the time domain waveform is measured and the anticausal part of the signal is artificially set to 0.

To handle properly the time domain wavefom in PyLayers, it is required to center the signal in the middle of the array.  The waveform has embedded in the object its frequency domain and time domain representation.

+ `st` member stands for signal in time domain
+ `sf` member stands for signal in frequency domain

```python
>>> print type(S.wav.sf)
>>> print type(S.wav.st)
```

+ `FUsignal` Frequency domain uniformly sampled base signal
+ `TUsignal` Time domain uniformly sampled base signal

## Construction of the propagation channel

The following representation shows the spatial spreading of the propagation channel.
On the left are scattered the intensity of rays wrt to angles of departure (in azimut and elevation).
On the right is the intensity of rays wrt to angles of arrival. It misses the application between the 2
planes as well as the delay dimension of the propagation channel.

```python
>>> from pylayers.antprop.signature import *
>>> from pylayers.antprop.channel import *
```

```python
>>> S.L.build()
```

```python
>>> S.L
```

```python
>>> S.L.Gt.pos
```

```python
>>> tx=np.array([759,1114,1.0])
>>> rx=np.array([767,1114,1.5])
>>> ctx = S.L.pt2cy(tx)
>>> crx = S.L.pt2cy(rx)
```

The sequence of command below :

+ initialize a signature between cycle ctx and cycle crx
+ evaluates the signature with a given cutoff value
+ calculates a set of 2D rays from signature and tx/rx coordinates
+ calculates a set of 3D ray from 2D rays and layout and ceil height (default H=3m)
+ calculates local basis and various geometric information out of the 3D ray and Layout
+ fill and reorganize the interactions object with proper material chararcteristics

```python
>>> Si = Signatures(S.L,ctx,crx)
>>> Si.run5(cutoff=5)
>>> r2d = Si.rays(tx,rx)
>>> r3d = r2d.to3D(S.L)
>>> r3d.locbas(S.L)
>>> r3d.fillinter(S.L)
```

Define a frequency base in GHz.

```python
>>> fGHz = np.arange(2,10,0.01)
```

Evaluate the propagation channel $\tilde{\mathbf{C}}$. Here the meaning of tilde is that the complex value of the channel do not include the phase term due to delay along the ray.

```python
>>> C = r3d.eval(fGHz)
```

## Construction of the transmission channel

The transmission channel is obtained from the combination of the propagation channel $\tilde{\mathbf{C}}$ and the vector antenna pattern at both side of the radio link. This operation is implemented in the `prop2tran` method of the `Ctilde` class.

```python
>>> sc = C.prop2tran()
```

The transmission channel is obtained by applying a vector radiation pattern using an antenna file.

In the presented case, it comes from a real antenna which has been used during the **FP7 project WHERE1** measurement campaign
M1.

```python
>>> sc
```

The antenna radiation pattern is stored in a very compact way thanks to Vector Spherical Harmonics decomposition.
The following gives information about the content of the antenna object.

```python
>>> S.tx.A.info()
```

The figure below plot on a same graph all the tansfer function in modulus and phase of the ray transfer function.

If a realistic antenna is applied it gives

```python
>>> sca = C.prop2tran(S.tx.A,S.rx.A)
```

## Calculate UWB Channel Impulse Response

Once the transmission channel has been evaluated on can convolved the waveform with the channel impulse response to get the received waveform.

```python
>>> r = sca.applywavB(S.wav.sfg)
```

```python
>>> r.y
```

```python
>>> fig,ax = r.plot(typ=['l20'])
>>> plt.axis([15,90,-120,-60])
>>> plt.title(u'Received Waveform $r(t)$')
```

```python
>>> r.plot(typ=['v'])
>>> #plt.axis([15,60,-0.3,0.3])
... plt.title(u'Received Waveform $r(t)$')
```

## Hermitian symetry enforcment

If the number of point for the transmission channel and the waveform were the same the mathematical operation is an Hadamrd-Shur product between
$\mathbf{Y}$ and $\mathbf{W}$.

$\mathbf{Y} = \mathbf{S} \odot \mathbf{W}$

In practice this is what is done after a resampling of the time base with a reinterpolated  time step.

The process which consists in going from time domain to frequency domain is delegated to a specialized class `pylayers.signal.bsignal.Bsignal` which maintains the proper
binding between signal samples and their indexation either in time or in frequency domain.

```python
>>> wgam = S.wav.sfg
>>> Y    = sc.apply(wgam)
>>> tau  = Y.taud
>>> dod = Y.dod
>>> doa = Y.doa
```

The transmission channel has a member data which is the time delay of each path in nano seconds. Notice that by default those delay are not sorted.

```python
>>> print 'tau =', tau[0:20]
```

```python
>>> h = plt.hist(tau,20)
```

Direction of arrival $(\theta_r,\phi_r)$ in radians

```python
>>> print "doa = ", doa[1:10,:]
```

```python
>>> plt.subplot(221)
>>> ht = plt.hist(doa[:,0],20)
>>> plt.xlabel(u'$\\theta_r$')
>>> plt.ylabel('#')
>>> plt.subplot(222)
>>> hp = plt.hist(doa[:,1],20)
>>> plt.xlabel(u'$\phi_r$')
>>> plt.ylabel('#')
>>> plt.subplot(223)
>>> ht = plt.hist(dod[:,0],20)
>>> plt.xlabel(u'$\\theta_t$')
>>> plt.ylabel('#')
>>> plt.subplot(224)
>>> hp = plt.hist(dod[:,1],20)
>>> plt.xlabel(u'$\phi_t$')
>>> plt.ylabel('#')
>>> plt.tight_layout()
```
