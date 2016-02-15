# Vector Spherical Harmonics Representation of Antennas


```python
>>> from pylayers.antprop.antenna import *
>>> from pylayers.antprop.antvsh import *
>>> %matplotlib inline
```

Loading an Antenna from a Matlab file

```python
>>> A = Antenna('S2R2.mat',directory='ant/UWBAN/Matfile')
```

The shape of the $F_{\phi}$ functions indicates :

- $N_f= 104$
- $N_{\theta} = 91$
- $N_{\phi} = 180 $

```python
>>> np.shape(A.Fp)
(91, 180, 104)
```

The frequency array is expressed in $GHz$ and delays are expressed in $ns$

```python
>>> fGHz = A.fGHz
```

```python
>>> fGHz.shape
(104,)
```

Then an electrical delay of $4.185ns$ is applied on the $F_{\theta}$

```python
>>> I = A.Ft[:,:,:]
```

```python
>>> I.shape
(91, 180, 104)
```

```python
...
>>> plt.figure(figsize=(10,8))
>>> plt.imshow(np.unwrap(np.angle(I[:,45,:])))
>>> plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
>>> plt.ylabel('f index')
>>> plt.colorbar()
>>> plt.figure()
>>> plt.plot(fGHz,np.unwrap(np.angle(I[45,85,:])))
>>> plt.xlabel('f index')
```

```python
>>> tau=4.185
>>> I = A.Ft[:,:,:]*np.exp(-2*1j*np.pi*fGHz[None,None,:]*tau)
```

```python
>>> plt.imshow(np.unwrap(np.angle(I[:,45,:])))
>>> plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
>>> plt.ylabel('f index')
>>> plt.colorbar()
>>> plt.figure()
>>> plt.plot(fGHz,np.unwrap(np.angle(I[45,85,:])))
[<matplotlib.lines.Line2D at 0x7f3f027bbe10>]
```

##### Display of the radiation pattern for all frequencies

```python
>>> plt.figure(figsize=(10,10))
>>> for nf in range(104):
...     plt.polar(A.phi,abs(A.Ft[45,:,nf]))
```

```python
>>> A.info()
 S2R2.mat
type :  mat
S2R2
Th1
04/13/12
09:59


2
2
Nb theta (lat) : 91
Nb phi (lon) : 180
No vsh coefficient calculated yet
```

# Evaluation of Vector Spherical Harmonics Coefficients

At that stage we compute the Vector Spherical Harmonics coefficients

```python
>>> A=vsh(A)
 
pass to vhaec
nlon =  180
nlat =  91
ityp =  0
nt =  1
idvw =  91
jdvw =  180
mdab =  90
ndab =  91
lvhaec =  40963
lwork =  57876
return from vhaec with br,bi,cr,ci
```

```python
>>> A.info()
S2R2.mat
type :  mat
S2R2
Th1
04/13/12
09:59


2
2
Nb theta (lat) : 91
Nb phi (lon) : 180
No vsh coefficient calculated yet
```

```python
>>> A.C.s1tos2(30)
```

```python
>>> A.C
```

```python
>>> fig = plt.figure(figsize=(8,8))
>>> A.C.show('s2',k=300)
```

```python
>>> A.C.s2tos3()
```

```python
>>> A.C
Br
-------------
L1  : 90
M1  : 89
Ncoeff s1 8010
NCoeff s2  : 495
Ncoeff s3 : 145

Bi
-------------
L1  : 90
M1  : 89
Ncoeff s1 8010
NCoeff s2  : 495
Ncoeff s3 : 145

Cr
-------------
L1  : 90
M1  : 89
Ncoeff s1 8010
NCoeff s2  : 495
Ncoeff s3 : 145

Ci
-------------
L1  : 90
M1  : 89
Ncoeff s1 8010
NCoeff s2  : 495
Ncoeff s3 : 145
```

```python
>>> fig = plt.figure(figsize=(8,8))
>>> A.C.show('s3')
>>> plt.tight_layout()
```

```python

```
