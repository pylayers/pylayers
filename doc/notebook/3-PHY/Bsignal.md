# Handling time and frequency domain signals : `Bsignal` Class

This section presents some features of the classes implemented in the [`pylayers.signal.bsignal.py`](http://pylayers.github.io/pylayers/modules/pylayers.signal.bsignal.html) module.

```python
>>> %matplotlib inline
```

The `Bsignal` class is a container for a signal with a base which can be either in time domain or frequency domain.

```python
>>> from pylayers.signal.bsignal import *
>>> from matplotlib.pyplot import *
```

As a first example, let construct an impulse signal normalized in energy. To do so there exist a specialized function : [`EnImpulse`](http://pylayers.github.io/pylayers/modules/generated/pylayers.signal.bsignal.EnImpulse.demo.html#pylayers.signal.bsignal.EnImpulse.demo)

```python
>>> E=EnImpulse(fe=40)
```

```python
>>> >>> from pylayers.signal.bsignal import *
>>> >>> ip    = EnImpulse(fc=4,band=3,thresh=10,fe=100)
>>> >>> Eip1  = ip.energy()
>>> >>> ESDu  = ip.esd(mode='unilateral')
>>> >>> ESDb  = ip.esd(mode='bilateral')
>>> >>> df    = ESDu.dx()
>>> >>> Eipu  = sum(ESDu.y)*df
>>> >>> Eipb  = sum(ESDb.y)*df
>>> >>> erru  = Eip1-Eipu
>>> >>> errb  = Eip1-Eipb
```

```python
>>> print Eip1
```

```python
>>> E.plot(typ='v')
```

```python
>>> E.energy()
```

The Fourier transform of this signal has the hermitian Symmetry.

```python
>>> F = E.fft()
>>> F.plot(typ='m')
```

```python
>>> F.y[0]
```

We then extract the non redundant part of the signal with the `ft` method

```python
>>> G=E.ft()
```

```python
>>> GH=G.symHz(100,scale='extract')
```

```python
>>> print GH.y[1]
>>> print GH.y[-1]
```

```python
>>> ip=F.ifft()
>>> ip2=GH.ifft()
```

```python
>>> f,a=E.plot(typ='v',labels=['original'])
>>> f,a=ip.plot(typ='v',fig=f,ax=a[0][0],labels=['no zero padding'])
>>> f,a=ip2.plot(typ='v',fig=f,ax=a[0][0],labels=['zero padding'])
>>> title('extract mode')
```

```python
>>> ip.energy()
```

```python
>>> ip2.energy()
```

```python
>>> Y=E.esd()
```



FHsignal for in CIR mode
------------------------

We create a Fusignal which corresponds to the signal

$$X_u(f) = \alpha e^{-2j\pi f \tau}$$

$$f\in [f_{min},f_{max}]$$

```python
>>> f = np.arange(2,10,0.01)
>>> y = 2*np.ones(len(f))*np.exp(-2*1j*np.pi*f*3)
>>> N = len(f)
>>> Hu = FUsignal(f,y)
>>> print N
```

```python
>>> Hu.plot(typ='m')
```

```python
>>> hu = Hu.ifft()
```

The inverse Fourier transform allows to recover perfectly the amplitude $\alpha$ and the delay $\tau$
of the channel

```python
>>> hu.plot(typ='m')
```

```python
>>> real=np.imag(hu.y)
>>> u = np.where(hu.y==max(hu.y))[0]
>>> tau = hu.x[u]
>>> alpha = abs(hu.y[u])
>>> print alpha,tau
```

```python
>>> H = Hu.symHz(100,scale='cir')
```

```python
>>> H.plot(typ='m')
```

```python
>>> h = H.ifft()
```

```python
>>> h.plot(typ='v')
```

```python
>>> real=np.imag(h.y)
>>> u = np.where(h.y==max(h.y))[0]
>>> tau = h.x[u]
>>> alpha = abs(h.y[u])
>>> print alpha,tau
```

```python
>>> fft.ifft(H.y)
```

```python
>>> print H.y[203]
>>> print H.y[-203]
>>> len(H.y)
```

```python
>>> Y=h.fft()
```

```python
>>> Y.plot(typ='m')
```
