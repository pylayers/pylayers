# Ctilde

```python
>>> %matplotlib inline
```

```python
>>> from pylayers.antprop.channel import *
>>> from pylayers.antprop.antenna import *
```

$\tilde{\mathbf{C}}$

```python
>>> C=Ctilde()
```

```python
>>> fGHz=np.r_[2.41,2.42]
```

```python
>>> pa = np.r_[0.018,-2.75,0.98]
>>> pb = np.r_[0.019,-2.75,0.98]
>>> C.los(pa=pa,pb=pb,fGHz=fGHz)
```

```python
>>> C
Ctilde
---------
(2,)
Nray : 1
fmin(GHz) : 2.41
fmax(GHz): 2.42
Nfreq : 2
```

```python
>>> print C.pa
>>> print C.pb
[ 0.018 -2.75   0.98 ]
[ 0.019 -2.75   0.98 ]
```

```python
>>> print C.tang
>>> print C.rang
[[ 1.57079633  0.        ]]
[[ 1.57079633 -3.14159265]]
```

# Antenna selection

```python
>>> Aa = Antenna('Omni')
>>> Ab = Antenna('3GPP_BackCenter_7.sh3')
```

```python
>>> Aa.polar()
(<matplotlib.figure.Figure at 0x7f3b71c65bd0>,
 <matplotlib.projections.polar.PolarAxes at 0x7f3b71dda8d0>)
```

```python
>>> Ab.polar(thd=90)
(<matplotlib.figure.Figure at 0x7f3b73c7cad0>,
 <matplotlib.projections.polar.PolarAxes at 0x7f3b71d07e50>)
```

```python
>>> Ta = np.eye(3)
>>> Tb = np.array([[-1,0,0],[0,-1,0],[0,0,1]])
```

```python
>>> Tb
array([[-1,  0,  0],
       [ 0, -1,  0],
       [ 0,  0,  1]])
```

```python
>>> C.tang.shape
(1, 2)
```

```python
>>> Cl=C.locbas(Tt=Ta,Tr=Tb)
```

```python
>>> H=Cl.prop2tran(a=Aa,b=Ab)
```

```python
>>> print H.Gap
>>> print H.Gat
>>> print H.Gbp
>>> print H.Gbt
[-inf]
[ 0.]
[-18.13467201+0.j]
[-7.86649093+0.j]
```

```python
>>> H.ak
array([ 3.9922806])
```

```python
>>> 3/0.3
10.0
```

```python
>>> C.tauk
array([ 0.00333333])
```

```python
>>> L=32.4+20*np.log10(2.4)+20*np.log10(d)
```

```python
>>> print L
[ 49.54664993]
```

```python
>>> ak=H.ak[0]
```

```python
>>> 20*np.log10(ak)
12.024421158892364
```

```python
>>> H
freq : 2.415 2.42 2
shape  : (1, 2)
tau (min, max) : 0.00333333333333 0.00333333333333
dist :0.001 0.001
Friis factor -j c/(4 pi f) has been applied
```

```python
>>> H.y
array([[  2.75228375e-15-3.99640698j,   2.74659722e-15-3.98814994j]])
```

```python
>>> ak
3.9922805962363754
```

```python
>>> H.Fat.y
array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
         1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]])
```

```python
>>> H.Fap.y
array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])
```

```python
>>> H.Fbt.y
array([[ 0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j,
         0.40427367 +2.78419052e-16j,  0.40427367 +2.78419052e-16j]])
```

```python
>>> C.Ctt.y
array([[1000.+0.j, 1000.+0.j]])
```

```python
>>> C.Ctp.y
array([[ 0.+0.j,  0.+0.j]])
```

```python
>>> H.alpha.y
array([[ 404.27366657 +2.78419052e-13j,  404.27366657 +2.78419052e-13j]])
```

```python
>>> H.y
array([[  2.75228375e-15-3.99640698j,   2.74659722e-15-3.98814994j]])
```

```python
>>> H.isFriis
True
```

```python
>>> 0.3/(4*np.pi*H.x)
array([ 0.0098854 ,  0.00986498])
```

```python

```

```python

```

```python

```
