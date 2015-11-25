---
!!python/unicode 'latex_envs':
  !!python/unicode 'bibliofile': !!python/unicode 'biblio.bib'
  !!python/unicode 'cite_by': !!python/unicode 'apalike'
  !!python/unicode 'current_citInitial': 1
  !!python/unicode 'eqLabelWithNumbers': true
  !!python/unicode 'eqNumInitial': 0
---

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
This is the first time the Layout is used. Graphs have to be built. Please Wait
```

```python
>>> L.show()
/home/uguen/anaconda/lib/python2.7/site-packages/matplotlib/collections.py:650: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors_original != str('face'):

/home/uguen/anaconda/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
(<matplotlib.figure.Figure at 0x7feaeb516290>,
 <matplotlib.axes._subplots.AxesSubplot at 0x7feaeb516a10>)
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
Antenna : Omni
Rotation matrice : 
 [[ 1.  0.  0.]
 [ 0.  1.  0.]
 [ 0.  0.  1.]]

Node b   
------  
position : [  761.0028967   1113.91576981     1.2       ]
Antenna : Omni
Rotation matrice : 
 [[ 1.  0.  0.]
 [ 0.  1.  0.]
 [ 0.  0.  1.]]

Link evaluation information : 
----------------------------- 
distance :  5.000 m 
delay : 16.667 ns
fmin (fGHz) : 2.4
fmax (fGHz) : 2.4
fstep (fGHz) : 0.0
Nf : 1
```

```python
>>> L.eval?
```

To evaluate a link there is the `eval` method. This method takes as argument 
+ a list of the desired outputs,
+ the type of algorithm being used, 
+ the ceil heigh 
+ the number of multi reflection between ceil and floor.

```python
>>> L.R.show(L=L.L,figsize=(10,10))
(<matplotlib.figure.Figure at 0x7feae5be9390>,
 <matplotlib.axes._subplots.AxesSubplot at 0x7feae5b44910>)
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
checkh5
Start Signatures
Signatures'> from 1_2_3 loaded
load signature
Stop signature 0.127700805664
Start Rays
Rays'> from 3_0_2 loaded
Stop rays 0.23652100563
Ctilde'> from 0_2_1 loaded
Tchannel'> from 0_2_1_0_0_0_0 loaded
```

```python
>>> plt.stem(aktk[1],aktk[0])
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
>>> f,a = C.show(cmap='jet',fig=fig,typ='l10',vmin=-100,vmax=-10)
/home/uguen/anaconda/lib/python2.7/site-packages/matplotlib/axes/_base.py:2562: UserWarning: Attempting to set identical left==right results
in singular transformations; automatically expanding.
left=2.4, right=2.4
  'left=%s, right=%s') % (left, right))
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
[<matplotlib.lines.Line2D at 0x7f17b5c56050>]
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
>>> dB=True
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
Trajectories performed in Layout : TA-Office.ini

Trajectory of agent John with ID 1
----------------------------------
t (s) : 0.00 : 0.20 : 119.80
dtot (m) : 443.20
Vmoy (m/s) : 3.70
                                 x         y        vx        vy        ax  \
t
1970-01-01 00:00:00.000  18.907750  2.528547  0.038749  0.155237  0.193744
1970-01-01 00:00:00.200  18.921699  2.584433  0.069748  0.279427  0.154995

                               ay      s
t
1970-01-01 00:00:00.000  0.776185  0.160
1970-01-01 00:00:00.200  0.620948  0.448

Trajectory of agent Alex with ID 2
----------------------------------
t (s) : 0.00 : 0.20 : 119.80
dtot (m) : 158.60
Vmoy (m/s) : 1.32
                                 x          y        vx        vy        ax  \
t
1970-01-01 00:00:00.000  24.306132  12.467593  0.030661 -0.157035  0.153303
1970-01-01 00:00:00.200  24.317170  12.411061  0.055189 -0.282663  0.122642

                               ay      s
t
1970-01-01 00:00:00.000 -0.785174  0.160
1970-01-01 00:00:00.200 -0.628139  0.448

Access point Router with ID 6
-----------------------------
t (s) : 0.00
Vmoy (m/s) : 0.0
              x  y    z  vx  vy  ax  ay  s
t
1970-01-01  0.5  2  2.5   0   0   0   0  0

Access point Router with ID 7
-----------------------------
t (s) : 0.00
Vmoy (m/s) : 0.0
              x   y    z  vx  vy  ax  ay  s
t
1970-01-01  0.7  14  2.5   0   0   0   0  0

Access point Router with ID 8
-----------------------------
t (s) : 0.00
Vmoy (m/s) : 0.0
             x   y    z  vx  vy  ax  ay  s
t
1970-01-01  39  13  2.5   0   0   0   0  0
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
