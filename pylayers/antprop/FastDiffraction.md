###This notebook is an implementation from the paper ["Fast UTD Diffraction Coefficient Using Only One Suare Root"](http://onlinelibrary.wiley.com/doi/10.1002/mop.28298/pdf) from Jean Franois Legendre and Thierry Marsault.

$N$ est compris entre 0 et 1. 0 est la lame de coupeau et 1 est le demi-plan

```python
>>> import numpy as np
>>> %pylab inline
>>> import  matplotlib.pyplot as plt
Populating the interactive namespace from numpy and matplotlib
```

```python
>>> fGHz = 2.4
>>> lambd = 0.3/fGHz
>>> k = 2*np.pi/lambd
>>> N = 0.5
>>> Nphi0 = 100
>>> Nphi = 50
>>> Nsi = 50
>>> Nsd = 60
```

$(\phi_0,\phi,s_i,s_d)$

```python
>>> si = np.linspace(0.1,4,Nsi)[np.newaxis,np.newaxis,:,np.newaxis]
>>> sd = np.linspace(0.1,10,Nsd)[np.newaxis,np.newaxis,np.newaxis,:]
>>> L = si*sd/(si+sd)
```

```python
>>> alpha = N*np.sqrt(L*2*k)
```

```python
>>> phi0 = np.linspace(0,N*np.pi,Nphi0)[:,np.newaxis,np.newaxis,np.newaxis]
>>> phi  = np.linspace(0,N*np.pi,Nphi)[np.newaxis,:,np.newaxis,np.newaxis]
```

```python
>>> xsi11 = (np.pi/4+phi-phi0)/(2*N)
>>> xsi12 = (np.pi/4-phi+phi0)/(2*N)
>>> xsi13 = (np.pi/4-phi-phi0)/(2*N)
>>> xsi14 = (phi+phi0-(2*N-1)*np.pi/4)/(2*N)
```

```python
>>> plt.figure(figsize=(10,10))
>>> subplot(141)
>>> imshow(xsi11[:,:,0,0])
>>> colorbar()
>>> subplot(142)
>>> imshow(xsi12[:,:,0,0])
>>> colorbar()
>>> subplot(143)
>>> imshow(xsi13[:,:,0,0])
>>> colorbar()
>>> subplot(144)
>>> imshow(xsi14[:,:,0,0])
>>> colorbar()
>>> tight_layout()
```

```python
>>> print xsi11.min()
>>> print xsi11.max()-np.pi/4
-0.785398163397
1.57079632679
```

```python
>>> print xsi12.min()
>>> print xsi12.max()-np.pi/4
-0.785398163397
1.57079632679
```

```python
>>> print xsi13.min()+np.pi/4
>>> print xsi13.max()
-1.57079632679
0.785398163397
```

```python
>>> print xsi14.min()
>>> print xsi14.max()-np.pi/4
0.0
2.35619449019
```

```python
>>> c11m = np.where(xsi11<0)
>>> c11p = np.where(xsi11>(np.pi/4))
>>> xsi11[c11m]=xsi11[c11m]*(-1)
>>> xsi11[c11p]=xsi11[c11p]-np.pi/4
```

```python
>>> c12m = np.where(xsi12<0)
>>> c12p = np.where(xsi12>(np.pi/4))
>>> xsi12[c12m]=xsi12[c12m]*(-1)
>>> xsi12[c12p]=xsi12[c12p]-np.pi/4
```

```python
>>> c13m = np.where(xsi13<0)
>>> xsi13[c13m]=xsi13[c13m]*(-1)
>>> c13p = np.where(xsi13>(np.pi/4))
>>> xsi13[c13p]=xsi13[c13p]-np.pi/4
```

```python
>>> #c14m = np.where(xsi14<0)
... #xsi14[c14m]=xsi14[c14m]*(-1)
... c14p = np.where(xsi14>(np.pi/4))
>>> xsi14[c14p]=xsi14[c14p]-np.pi/4
```

```python
>>> plt.subplot(141)
>>> imshow(xsi11[:,:,0,0])
>>> colorbar()
>>> plt.subplot(142)
>>> imshow(xsi12[:,:,0,0])
>>> colorbar()
>>> plt.subplot(143)
>>> imshow(xsi13[:,:,0,0])
>>> colorbar()
>>> plt.subplot(144)
>>> plt.imshow(xsi14[:,:,0,0])
>>> colorbar()
```

```python
>>> np.where(xsi14>np.pi/4)
(array([ 1,  2,  3, ..., 99, 99, 99]),
 array([49, 49, 48, ..., 47, 48, 49]),
 array([0, 0, 0, ..., 0, 0, 0]),
 array([0, 0, 0, ..., 0, 0, 0]))
```

```python
>>> y11 = alpha * np.sin(xsi11)
>>> y12 = alpha * np.sin(xsi12)
>>> y13 = alpha * np.sin(xsi13)
>>> y14 = alpha * np.sin(xsi14)
```

```python
>>> st11 = np.ones(np.shape(xsi11)[0:2])[:,:,np.newaxis,np.newaxis]
>>> st11[c11m]=-st11[c11m]
>>> st11[c11p]=-st11[c11p]
>>> st12 = np.ones(np.shape(xsi12)[0:2])[:,:,np.newaxis,np.newaxis]
>>> st12[c12m]=-st12[c12m]
>>> st12[c12p]=-st12[c12p]
>>> st13 = np.ones(np.shape(xsi13)[0:2])[:,:,np.newaxis,np.newaxis]
>>> st13[c13m]=-st13[c13m]
>>> st13[c13p]=-st13[c13p]
>>> st14 = np.ones(np.shape(xsi14)[0:2])[:,:,np.newaxis,np.newaxis]
>>> #st14[c14m]=-st14[c14m]
... st14[c14p]=-st14[c14p]
```

```python
>>> subplot(141)
>>> imshow(st11[:,:,0,0])
>>> subplot(142)
>>> imshow(st12[:,:,0,0])
>>> subplot(143)
>>> imshow(st13[:,:,0,0])
>>> subplot(144)
>>> imshow(st14[:,:,0,0])
```

```python
>>> np.where(y11*y11<0.005)
(array([ 0,  0,  0, ..., 99, 99, 99]),
 array([ 1,  1,  1, ..., 26, 26, 26]),
 array([0, 0, 0, ..., 4, 5, 6]),
 array([0, 1, 2, ..., 0, 0, 0]))
```

```python
>>> y11c = y11*y11
>>> cinf = np.where(y11c<0.005)
>>> csup = np.where(y11c>=0.005)
>>> yinf = y11[cinf]
>>> ysup = y11[csup]
>>> ginf = -0.5+0.39384228*yinf*(1+1j)
>>> #tmp = 0.5/(yp*yp)
... tmp = 0.5*y11c[csup]
>>> gsup = -(0.1994711402/ysup)*(tmp+1j*(tmp-1))
```

```python
>>> np.shape(y11)
(100, 50, 50, 60)
```

```python
>>> np.shape(ginf)
(296711,)
```

```python
>>> np.shape(gsup)
(14703289,)
```

```python
>>> 100*50*50*60 - (14703289+296711)
0
```

```python
>>> g11 = np.zeros(np.shape(y11),dtype=complex)
```

```python
>>> np.shape(g11)
(100, 50, 50, 60)
```

```python
>>> np.shape(cinf)
(4, 296711)
```

```python
>>> g11[cinf]=ginf
>>> g11[csup]=gsup
```

```python
>>> #g11
```

```python
>>> def g(y):
...     if y*y<0.5:
...         g = -0.5+0.39384228*y*(1+1j)
...     else:
...         g = (-0.1994711402/y)*((0.5/y*y)+1)+1j*((0.5/y*y) - 1)
...     return (g)
```

```python
>>> y = np.arange(-10,10,0.001)
```

```python
>>> f = []
>>> for x in y:
...     f.append(g(x))
```

```python
>>> plt.plot(y,np.abs(np.array(f)))
[<matplotlib.lines.Line2D at 0x7ff2d0ada750>]
```

```python
>>> from simpy  import *
```

```python
>>> 1/0.033333
30.000300003000028
```

```python
>>> def sinq21(xsi21):
...     s = xsi21*((1-1/6)*xsi21*xsi21*(1-0.05*xsi21*xsi21*(1-0.1428571429*xsi21*xsi21)))
...     return s
```

```python
>>> s
<function __main__.sinq21>
```

```python
>>> a
<function __main__.sinq21>
```

```python
>>> 1/0.0178571428
56.0000001792
```

```python
>>> 1/0.05
20.0
```

```python

```
