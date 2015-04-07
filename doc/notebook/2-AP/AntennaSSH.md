# Scalar Spherical Harmonics

```python
>>> from pylayers.antprop.antenna import *
>>> from pylayers.antprop.antssh import *
>>> %matplotlib inline
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
```

```python
>>> A = Antenna('S1R1.mat',directory='ant/UWBAN/Matfile')
```

```python
>>> A
FileName : S1R1.mat
-----------------------
fmin : 0.80GHz
fmax : 5.95GHz
step : 50.00MHz
Nf : 104
-----------------------
Ntheta : 91
Nphi : 180
GmaxdB : 2.23 dB
   f = 5.60 GHz
   theta = 70.00 (degrees)
   phi = 272.00  (degrees)
antenna name : Th1
date : 04/12/12
time : 15:55
Notes : Mohamed at the log

Serie : 1
Run : 1
Nb theta (lat) : 91
Nb phi (lon) :180
```

To calculate scalar spherical harmonics use method `ssh(A,L)`

```python
>>> L = 5
>>> A = ssh(A,L=5)
```

```python
>>> A
FileName : S1R1.mat
-----------------------
fmin : 0.80GHz
fmax : 5.95GHz
step : 50.00MHz
Nf : 104
-----------------------
Ntheta : 91
Nphi : 180
GmaxdB : 2.23 dB
   f = 5.60 GHz
   theta = 70.00 (degrees)
   phi = 272.00  (degrees)
antenna name : Th1
date : 04/12/12
time : 15:55
Notes : Mohamed at the log

Serie : 1
Run : 1
Nb theta (lat) : 91
Nb phi (lon) :180
```

```python
>>> plt.plot(abs(A.S.Cx.s2[0]))
[<matplotlib.lines.Line2D at 0x7f7e259c47d0>]
```

```python
>>> A.savesh2()
/home/uguen/Bureau/P1/ant/S1R1.sh2  already exist
```

```python
>>> A.loadsh2()
```

```python
>>> plt.plot(abs(A.S.Cx.s2[0]))
[<matplotlib.lines.Line2D at 0x7f7e258c67d0>]
```

```python
>>> A.S.s2tos3()
```

```python
>>> plt.plot(abs(A.S.Cx.s3[0]))
[<matplotlib.lines.Line2D at 0x7f7e25815750>]
```

```python
>>> A.S.Cx.ind2.shape
(36, 2)
```

```python
>>> A.savesh3()
/home/uguen/Bureau/P1/ant/S1R1.sh3  already exist
```

```python
>>> plt.plot(abs(A.S.Cx.s2[0]))
[<matplotlib.lines.Line2D at 0x7f7e25752d90>]
```

```python
>>> A.loadsh3()
```

```python
>>> plt.plot(abs(A.S.Cx.s3[100]))
[<matplotlib.lines.Line2D at 0x7f7e2569e410>]
```

```python
>>> plt.plot(abs(A.S.Cx.s2[100]))
[<matplotlib.lines.Line2D at 0x7f7e25569290>]
```

```python
>>> A.__dict__.keys()
['tau',
 'Nf',
 'PhotoFile',
 'Np',
 'Nt',
 'Run',
 'source',
 '_filename',
 'Serie',
 'Ftheta',
 'theta',
 'fromfile',
 'phi',
 'Fphi',
 'Notes',
 'fa',
 'S',
 'AntennaName',
 'typ',
 'DataFile',
 'evaluated',
 'Date',
 'ext',
 'StartTime',
 'SqG']
```

```python
>>> A.S.Cx.__dict__.keys()
['k2', 'ind3', 'ind2', 'fmax', 's2', 'Nf', 's3', 'lmax', 'fmin']
```

```python
>>> A.S.Cx
Nf   : 104
fmin (GHz) : 0.8
fmax (GHz) : 5.95
NCoeff s2  : 36
Ncoeff s3 : 143
```
