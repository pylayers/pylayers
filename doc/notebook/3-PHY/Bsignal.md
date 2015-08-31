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
(<matplotlib.figure.Figure at 0x7fd3bcafb690>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fd3bcaf3c50>]], dtype=object))
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
(<matplotlib.figure.Figure at 0x7fd3bc812a10>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fd3bc8254d0>]], dtype=object))
```

```python
>>> h = H.ifft()
```

```python
>>> h.plot(typ='v')
(<matplotlib.figure.Figure at 0x7fd3bc972790>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fd3bcaa1310>]], dtype=object))
```

```python
>>> real=np.imag(h.y)
>>> u = np.where(h.y==max(h.y))[0]
>>> tau = h.x[u]
>>> alpha = abs(h.y[u])
>>> print alpha,tau
[[  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
    1.91361122e-02   1.16994371e-03]
 [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
    1.91361122e-02   1.16994371e-03]
 [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
    1.91361122e-02   1.16994371e-03]
 ..., 
 [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
    1.91361122e-02   1.16994371e-03]
 [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
    1.91361122e-02   1.16994371e-03]
 [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
    1.91361122e-02   1.16994371e-03]] [-49.97728305 -49.97728305 -49.97728305 ..., -49.97728305 -49.97728305
 -49.97728305]
```

```python
>>> fft.ifft(H.y)
array([[ -1.93565190e-15 -1.70240923e-19j,
          2.62295322e-04 -3.27871407e-19j,
          8.73458329e-04 -4.09839258e-20j, ...,
         -1.06670199e-04 +2.90350482e-19j,
         -8.69428086e-04 -1.58117458e-18j,
         -5.31550980e-05 -2.71727936e-20j]])
```

```python
>>> print H.y[...,203]
>>> print H.y[...,-203]
>>> len(H.y)
[ 0.11616926-0.04599466j]
[ 0.11616926+0.04599466j]
1
```

```python
>>> Y=h.fft()
```

```python
>>> Y.plot(typ='m')
(<matplotlib.figure.Figure at 0x7fd3bc8bd8d0>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fd3bcd15dd0>]], dtype=object))
```
