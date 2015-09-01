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
>>> E=TUsignal()
>>> E.EnImpulse(feGHz=40)
```

```python
>>> E.plot(typ='v')
(<matplotlib.figure.Figure at 0x7f8b08493810>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad8a32b10>]], dtype=object))
```

```python
>>> E.energy()
array([ 1.00000008])
```

The Fourier transform of this signal exhibits the Hermitian Symmetry.

```python
>>> F = E.fft()
>>> F.plot(typ='m')
(<matplotlib.figure.Figure at 0x7f8ad89e0c50>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad8a05e90>]], dtype=object))
```

We then extract the non redundant part of the signal with the `ft` method

```python
>>> G=E.ft()
```

```python
>>> GH=G.symHz(100,scale='extract')
```

```python
>>> print GH.y[0,1]
>>> print GH.y[0,-1]
(-0.0014441784194-4.88037298122e-05j)
(-0.0014441784194+4.88037298122e-05j)
```

```python
>>> ip = F.ifft()
>>> ip2= GH.ifft()
```

```python
>>> f,a=E.plot(typ='v',labels=['original'])
>>> f,a=ip.plot(typ='v',fig=f,ax=a[0][0],labels=['no zero padding'])
>>> f,a=ip2.plot(typ='v',fig=f,ax=a[0][0],labels=['zero padding'])
>>> title('extract mode')
```

```python
>>> ip.energy()
array([ 1.00000008])
```

```python
>>> ip2.energy()
array([ 3.18478273])
```

```python
>>> Y=E.esd()
```

FHsignal in CIR mode
------------------------

We create a Fusignal which corresponds to the signal

$$X_u(f) = \alpha e^{-2j\pi f \tau}$$

$$f\in [f_{min},f_{max}]$$

```python
>>> fGHz = np.arange(2,10,0.01)
>>> tau = 20
>>> y = 2*np.ones(len(fGHz))*np.exp(-2*1j*np.pi*fGHz*tau)
>>> Hu = FUsignal(fGHz,y)
```

```python
>>> Hu.plot(typ='m')
>>> Hu.plot(typ='r')
(<matplotlib.figure.Figure at 0x7f8ad84b6050>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad84630d0>]], dtype=object))
```

```python
>>> hu = Hu.ifft()
```

The inverse Fourier transform allows to recover perfectly the amplitude $\alpha$ and the delay $\tau$
of the channel

```python
>>> hu.plot(typ='m')
(<matplotlib.figure.Figure at 0x7f8ad838ec90>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad84b6c90>]], dtype=object))
```

```python
>>> real=np.imag(hu.y)
>>> u = np.where(hu.y==max(hu.y))[0]
>>> tau = hu.x[u]
>>> alpha = abs(hu.y[u])
```

```python
>>> H = Hu.symHz(100,scale='cir')
```

```python
>>> H.plot(typ='m')
(<matplotlib.figure.Figure at 0x7f8ad8282450>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad8312f10>]], dtype=object))
```

```python
>>> h = H.ifft()
```

```python
>>> h.plot(typ='v')
(<matplotlib.figure.Figure at 0x7f8ad81c6890>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad83565d0>]], dtype=object))
```

```python
>>> real=np.imag(h.y)
>>> u = np.where(h.y==max(h.y))[0]
>>> tau = h.x[u]
>>> alpha = abs(h.y[u])
```

```python
>>> fft.ifft(H.y)
array([[ -1.50593859e-15 -6.41964563e-20j,
          1.22745263e-04 -1.36427337e-19j,
          8.94216494e-05 -1.03247967e-19j, ...,
          1.05839739e-05 +7.80645228e-20j,
         -1.37135712e-04 -1.94405223e-19j,
          8.17123442e-05 +3.02799103e-19j]])
```

```python
>>> print H.y[...,203]
>>> print H.y[...,-203]
>>> len(H.y)
[-0.10108118-0.07343977j]
[-0.10108118+0.07343977j]
1
```

```python
>>> Y=h.fft()
```

```python
>>> Y.plot(typ='m')
(<matplotlib.figure.Figure at 0x7f8ad836ec90>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f8ad8491590>]], dtype=object))
```
