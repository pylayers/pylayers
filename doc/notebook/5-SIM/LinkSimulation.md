```python
>>> from pylayers.simul.link import *
>>> %matplotlib inline
>>> import seaborn as sns
>>> sns.set_style("dark")
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
```

# How to do Ray Tracing simulation using DLink

This section illustrates the [`link`](http://pylayers.github.io/pylayers/modules/pylayers.simul.link.html)  module.
A `Dlink` object is a deterministic (Single Input Single Output) SISO link.

```python
>>> L=DLink()
```

```python
>>> L.show()
(<matplotlib.figure.Figure at 0x7f7b8d6aa910>,
 <matplotlib.axes._subplots.AxesSubplot at 0x7f7b8d6c3a50>)
```

```python
>>> L
filename: Links_0_defstr.ini.h5
Link Parameters :
------- --------
Layout : defstr.ini

Node a   
------  
position : [  766.00300113  1113.94747911     1.2       ]
Antenna : S2R2.sh3
Rotation matrice : 
 [[ 1.  0.  0.]
 [ 0.  1.  0.]
 [ 0.  0.  1.]]

Node b   
------  
position : [  761.0028967   1113.91576981     1.2       ]
Antenna : S2R2.sh3
Rotation matrice : 
 [[ 1.  0.  0.]
 [ 0.  1.  0.]
 [ 0.  0.  1.]]

Link evaluation information : 
----------------------------- 
distance :  5.000 m 
delay : 16.667 ns
fmin (fGHz) : 2.0
fmax (fGHz) : 11.0
fstep (fGHz) : 0.05
Nf : 181
```

```python
>>> L.eval()
Signatures'> from 1_2_3 loaded
Rays'> from 3_0_2 loaded
Ctilde'> from 0_2_0 loaded
Tchannel'> from 0_2_0_0_0_1_1 loaded
(array([  8.74446418e-05,   1.46698690e-04,   0.00000000e+00,
          3.14995445e-06,   1.51702810e-04,   7.29356204e-05,
          5.90710148e-05,   2.46454264e-05,   2.44450202e-05,
          4.87319705e-06,   1.75686332e-04,   7.13675689e-06,
          7.47211694e-05,   6.74741489e-08,   1.95936070e-06,
          2.23741852e-07,   2.26646476e-06,   4.14435099e-06,
          3.00650318e-06,   3.00698950e-06,   2.61297158e-05,
          2.70506446e-05,   0.00000000e+00,   0.00000000e+00,
          7.97779088e-08,   2.38562292e-06,   4.04156289e-08,
          2.64709912e-06,   9.07630898e-08,   1.12916970e-06,
          1.22492177e-07,   1.33692589e-06,   8.04381204e-06,
          8.49458480e-06,   6.14904213e-05,   5.66523749e-05,
          4.67849894e-06,   3.68188625e-06,   3.68423370e-06,
          0.00000000e+00,   3.88676042e-06,   3.88632416e-06,
          2.23404583e-07,   5.46676503e-06,   1.02274728e-05,
          1.00586800e-05,   1.07393894e-05,   9.71108468e-06,
          9.25607532e-08,   1.62695253e-07,   4.71935354e-07,
          1.21853846e-06,   1.61900680e-07,   1.01020395e-07,
          1.47301403e-06,   6.52210803e-07,   4.42388951e-06,
          6.41935326e-06,   5.64142360e-06,   3.84647619e-06,
          3.84739961e-06,   5.64145807e-06,   3.27350264e-06,
          3.34549293e-06,   9.67273293e-08,   1.64206449e-07,
          2.35671508e-07,   6.21878652e-06,   1.36855910e-07,
          3.45202976e-06,   2.66614869e-07,   9.25499264e-07,
          9.77377034e-07,   0.00000000e+00,   0.00000000e+00,
          1.01308828e-07,   1.36303838e-07,   0.00000000e+00,
          0.00000000e+00,   1.91570108e-07,   2.77178941e-07,
          3.40906186e-06,   5.08162922e-06,   2.75700486e-07,
          0.00000000e+00,   1.71902246e-06,   1.71676236e-06,
          1.68408980e-06,   1.76350549e-06,   1.12834353e-07,
          3.06875167e-08,   1.41120249e-07,   3.95694440e-08,
          3.54499281e-07,   4.25309327e-07]),
 array([ 16.66734994,  18.48784882,  20.53778357,  33.33365909,
         33.3343511 ,  26.03460301,  26.03460301,  31.44237961,
         31.44261404,  34.28021045,  34.28088335,  35.42785385,
         35.42850495,  35.89995802,  35.89997275,  35.90181477,
         35.90184832,  50.00066032,  50.0013426 ,  50.00136179,
         32.44415565,  32.44438284,  33.65446829,  33.65468731,
         36.78052454,  36.78053892,  36.78233685,  36.78236959,
         37.85243698,  37.85245096,  37.85419797,  37.85422979,
         38.87329197,  38.87329197,  38.87388536,  38.87388536,
         50.63660763,  50.63728134,  50.63730029,  51.42048261,
         51.42114606,  51.42116471,  66.66765177,  66.66836297,
         37.26423534,  37.26423534,  37.26443315,  37.26443315,
         41.09509686,  41.09509686,  41.09510973,  41.09510973,
         41.0967189 ,  41.0967189 ,  41.0967482 ,  41.0967482 ,
         53.85226116,  53.85226116,  53.85289465,  53.85289465,
         53.85291246,  53.85291246,  55.86686011,  55.8671227 ,
         56.66325336,  56.66785273,  67.14592909,  67.14663522,
         67.7390271 ,  67.73972705,  83.33535519,  56.43674387,
         56.43700381,  57.14110656,  57.14136329,  57.2252067 ,
         57.2297609 ,  57.91998171,  57.92448129,  69.60298696,
         69.60298696,  69.60366816,  69.60366816,  83.71846526,
         84.19490142,  59.33890847,  59.33890847,  59.3391557 ,
         59.3391557 ,  60.08930256,  60.08930256,  60.09363971,
         60.09363971,  85.70170024,  85.70170024]))
```

To evaluate a link there is the `eval` method. This method takes as argument 
+ a list of the desired outputs,
+ the type of algorithm being used, 
+ the ceil heigh 
+ the number of multi reflection between ceil and floor.

```python
>>> L.R.show3()
```

```python
>>> L.H.taud
array([ 16.66734994,  18.48784882,  20.53778357,  33.33365909,
        33.3343511 ,  26.03460301,  26.03460301,  31.44237961,
        31.44261404,  34.28021045,  34.28088335,  35.42785385,
        35.42850495,  35.89995802,  35.89997275,  35.90181477,
        35.90184832,  50.00066032,  50.0013426 ,  50.00136179,
        32.44415565,  32.44438284,  33.65446829,  33.65468731,
        36.78052454,  36.78053892,  36.78233685,  36.78236959,
        37.85243698,  37.85245096,  37.85419797,  37.85422979,
        38.87329197,  38.87329197,  38.87388536,  38.87388536,
        50.63660763,  50.63728134,  50.63730029,  51.42048261,
        51.42114606,  51.42116471,  66.66765177,  66.66836297,
        37.26423534,  37.26423534,  37.26443315,  37.26443315,
        41.09509686,  41.09509686,  41.09510973,  41.09510973,
        41.0967189 ,  41.0967189 ,  41.0967482 ,  41.0967482 ,
        53.85226116,  53.85226116,  53.85289465,  53.85289465,
        53.85291246,  53.85291246,  55.86686011,  55.8671227 ,
        56.66325336,  56.66785273,  67.14592909,  67.14663522,
        67.7390271 ,  67.73972705,  83.33535519,  56.43674387,
        56.43700381,  57.14110656,  57.14136329,  57.2252067 ,
        57.2297609 ,  57.91998171,  57.92448129,  69.60298696,
        69.60298696,  69.60366816,  69.60366816,  83.71846526,
        84.19490142,  59.33890847,  59.33890847,  59.3391557 ,
        59.3391557 ,  60.08930256,  60.08930256,  60.09363971,
        60.09363971,  85.70170024,  85.70170024])
```

```python
>>> aktk=L.eval(force=[], output=['sig','ray','Ct','H'],
...             si_algo='old',ra_ceil_height_meter=3,ra_number_mirror_cf=1)
Signatures'> from 1_2_3 loaded
Rays'> from 3_0_2 loaded
Ctilde'> from 0_2_0 loaded
Tchannel'> from 0_2_0_0_0_1_1 loaded
```

```python
>>> plt.stem(aktk[1],aktk[0])
<Container object of 3 artists>
```

The propagation channel (without antenna) can be vizualized on a ray by ray mode.

```python
>>> type(L.C)
pylayers.antprop.channel.Ctilde
```

```python
>>> #L._show3()sns.set_style("dark")
```

```python
>>> fig = plt.figure(figsize=(8,8))
>>> C = L.C
>>> f,a = C.show(cmap='jet',fig=fig)
```

It is possible to look at individual ray transfer function, as illustrated below.

```python
>>> C.Ctt.y.shape
(95, 181)
```

```python
>>> ir = 80
>>> plt.plot(C.Ctt.x,abs(C.Ctt.y[ir,:]))
>>> plt.xlabel('Frequency (GHz)')
>>> plt.ylabel('Level (linear)')
>>> plt.title('Modulus of the ray '+str(ir)+' transfer function')
```

```python
>>> ir = 30
>>> plt.plot(C.Ctt.x,abs(C.Ctt.y[ir,:]))
>>> plt.xlabel('Frequency (GHz)')
>>> plt.ylabel('Level (linear)')
>>> plt.title('Modulus of the ray '+str(ir)+' transfer function')
```

In the link we also have the transmission channel accounting for the effect of antennas and Friis factor. If the ray transfer function is scaled with $\frac{4\pi f}{c}$

```python
>>> plt.plot(L.H.x,L.H.y[0,:]*4*np.pi*L.H.x/0.3)
/home/uguen/anaconda/lib/python2.7/site-packages/numpy/core/numeric.py:462: ComplexWarning: Casting complex values to real discards the imaginary part
  return array(a, dtype, copy=False, order=order)
[<matplotlib.lines.Line2D at 0x7f7b8f877090>]
```

Notice that in this case the frequency

The infinite bandwidth channel impulse response is plotted below from the extrated set $\{\alpha_k,\tau_k\}$.

```python
>>> plt.stem(aktk[1],aktk[0])
>>> plt.title('Infinite bandwith Channel Impulse response')
>>> plt.xlabel('delay (ns)')
>>> plt.ylabel('amplitude (linear scale')
```

```python
>>> import pylayers.simul.simulnet as sn
>>> import pylayers.simul.simultraj as st
```

```python
>>> S=sn.Simul()
Layout graphs are loaded from /home/uguen/Bureau/P1/struc/ini
```

```python
>>> S.L

----------------
TA-Office.ini
Image('/home/uguen/Bureau/P1/struc/images/DLR4991.png')
----------------

Number of points  : 71
Number of segments  : 87
Number of sub segments  : 16
Number of cycles  : 18
Number of rooms  : 17
degree 0 : []
degree 1 : []
number of node point of degree 2 : 39
number of node point of degree 3 : 32

xrange :(0.0, 40.0)
yrange :(0.0, 15.0)

Useful dictionnaries
----------------
dca {cycle : []} cycle with an airwall
sl {slab name : slab dictionary}
name :  {slab :seglist} 

Useful arrays
----------------
pt : numpy array of points 
normal : numpy array of normal 
offset : numpy array of offset 
tsg : get segment index in Gs from tahe
isss :  sub-segment index above Nsmax
tgs : get segment index in tahe from Gs
lsss : list of segments with sub-segment
sla : list of all slab names (Nsmax+Nss+1)
degree : degree of nodes
```

```python
>>> S.runsimul()
```

```python
>>> from pylayers.mobility.trajectory import *
```

A trajectories object is a list of trajectories. The loadh5 methods by default loads the file '.h5'
generated by [`Simulnet`](http://pylayers.github.io/pylayers/modules/pylayers.simul.simulnet.html).

```python
>>> T=Trajectories()
>>> T.loadh5()
```

```python
>>> T
```

A SimulTraj object is derived from a trajectory calculated previously in simulnet and a body agent description.
The Simultraj object get the trajectories from the `simultaj.ini` file.

```python
>>> St=st.Simul(verbose=False)
```

```python
>>> #St.run(t=list(np.arange(0,1,0.1)),OB=True,B2B=True,B2I=True)
```

```python
>>> #St.data
```

Information about the simulated network is obtained

```python
>>> St.N
```

```python
>>> #St._show3()
```

```python
>>> #St.data.head()
```

```python
>>> #ak,tk,ek=St._loadh5(2,'0_Alex','1_Alex','bluetooth-class2')
```

```python
>>> #stem(tk,ak)
```

```python

```
