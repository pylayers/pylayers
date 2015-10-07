# Effect of Modyfiying the Nature of Sub-Segments

This notebook illustrtates a simple ray tracing simulation with diffecent material properties for a single segment separating 2 rooms which contains multi-subsegments. The noteboook illustrates in details the whole steps.

```python
>>> from pylayers.simul.link import *
>>> from pylayers.antprop.rays import *
>>> from pylayers.gis.layout import *
>>> from pylayers.antprop.signature import *
>>> import pylayers.signal.bsignal as bs
>>> import pylayers.signal.waveform as wvf
>>> from pylayers.simul.simulem import *
>>> import matplotlib.pyplot as plt
>>> %matplotlib inline
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
Warning : OSM Parser seems to be not installed
```

Let start by loading a simple layout with 2 single rooms. The multi subsegment appears in the middle with the red vertical lines. Each subsegment is materialized by a  segment.

```python
>>> L=Layout('defstr.ini')
>>> f,a=L.showG('s',subseg=True,figsize=(10,10))
```

The studied configuration is composed of a simple 2 rooms building separated by a subsegment which has a multi subsegment attribute. The attribute of the subsegment can be changed  with the method [`chgmss`](http://pylayers.github.io/pylayers/modules/generated/pylayers.gis.layout.Layout.chgmss.html) (change multisubsegment). In the example WOOD in the lower part then 10cm of AIR then wood again until the ceiling.

```python
>>> L.chgmss(1,ss_name=['WOOD','AIR','WOOD'],ss_z =[(0.0,2.7),(2.7,2.8),(2.8,3)],ss_offset=[0,0,0])
```

As the Layout structure has been modified it is required to rebuild the structure.

```python
>>> L.build()
>>> L.save()
structure saved in  defstr.str2
structure saved in  defstr.ini
```

The $\mathcal{G}_s$ graph dictionnary has the following structure

```python
>>> L.Gs.node
{-8: {},
 -7: {},
 -6: {},
 -5: {},
 -4: {},
 -3: {},
 -2: {},
 -1: {},
 1: {'connect': [-8, -7],
  'name': 'PARTITION',
  'ncycles': [1, 2],
  'norm': array([-0.999982  , -0.00599989,  0.        ]),
  'offset': 0,
  'ss_name': ['WOOD', 'AIR', 'WOOD'],
  'ss_offset': [0, 0, 0],
  'ss_z': [(0.0, 2.7), (2.7, 2.8), (2.8, 3)],
  'transition': True,
  'z': (0.0, 3.0)},
 2: {'connect': [-8, -2],
  'name': 'WALL',
  'ncycles': [1, 2],
  'norm': array([ 0.99997778,  0.00666652,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 3: {'connect': [-7, -5],
  'name': 'WALL',
  'ncycles': [1, 2],
  'norm': array([-0.99997775, -0.00667097,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 4: {'connect': [-6, -1],
  'name': 'WALL',
  'ncycles': [2, 0],
  'norm': array([ 0.99997888,  0.00649986,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 5: {'connect': [-6, -5],
  'name': 'WALL',
  'ncycles': [2, 0],
  'norm': array([-0.00619988,  0.99998078,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 6: {'connect': [-5, -4],
  'name': 'WALL',
  'ncycles': [1, 0],
  'norm': array([-0.00639987,  0.99997952,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 7: {'connect': [-4, -3],
  'name': 'WALL',
  'ncycles': [1, 0],
  'norm': array([ 0.99997887,  0.00650149,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 8: {'connect': [-3, -2],
  'name': 'WALL',
  'ncycles': [1, 0],
  'norm': array([ 0.00639987, -0.99997952,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)},
 9: {'connect': [-2, -1],
  'name': 'WALL',
  'ncycles': [2, 0],
  'norm': array([ 0.00639987, -0.99997952,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)}}
```

We define now two points which are the termination of a radio link.

```python
>>> #tx=np.array([759,1114,1.5])
... #rx=np.array([767,1114,1.5])
... tx=np.array([759,1114,1.5])
>>> rx=np.array([767,1114,1.5])
```

```python
>>> L.chgmss(1,ss_name=['WOOD','AIR','WOOD'],ss_z =[(0.0,2.7),(2.7,2.8),(2.8,3)],ss_offset=[0,0,0])
>>> L.save()
>>> fGHz=np.linspace(1,11,100)
>>> #Aa = Antenna('S1R1.vsh3')
... #Ab = Antenna('S1R1.vsh3')
... Aa = Antenna('Gauss',fGHz=fGHz)
>>> Ab = Antenna('Gauss',fGHz=fGHz)
>>> Lk = DLink(L=L,a=tx,b=rx,Aa=Aa,Ab=Ab,fGHz=np.linspace(1,11,100))
structure saved in  defstr.str2
structure saved in  defstr.ini
```

A link is the set of a layout and 2 termination points.

```python
>>> f,a=Lk.show(rays=True)
```

On the figure above, we can see the Tx and Rx each placed in a different room appart from a wall with a subsegement placed in the middle.
Then for evaluating the radio link, simply type:

```python
>>> ak,tauk=Lk.eval(force=True,a=tx,b=rx,applywav=True)
checkh5
Start Signatures
algo 7
Signatures'> from 2_1_3 saved
Stop signature 0.0473918914795
Start Rays
Rays'> from 3_2_3 saved
Stop rays 0.465720891953
Ctilde'> from 2_3_0 saved
```

```python
>>> Lk.C
Ctilde
---------
(155, 100)
Nray : 155
fmin(GHz) : 1.0
fmax(GHz): 11.0
Nfreq : 100
```

```python
>>> f = plt.figure(figsize=(10,10))
>>> f,a=Lk.C.show(cmap='jet',fig=f,typ='l20',vmin=-120,vmax=-10)
```

```python
>>> fGHz=np.arange(2,6,0.5)
>>> wav = wvf.Waveform(fcGHz=4,bandGHz=1.5)
>>> wav.show()
```

```python
>>> wav.st.y.shape
(1, 251)
```

```python
>>> len(Lk.fGHz)
100
```

```python
>>> Lk = DLink(L=L,a=tx,b=rx)
```

```python
>>> Lk.a
array([  759. ,  1114. ,     1.5])
```

```python
>>> Lk.b
array([  767. ,  1114. ,     1.5])
```

```python
>>> cir = Lk.H.applywavB(wav.sf)
```

```python
>>> layer = ['AIR','AIR','AIR']
>>> Lk.L.chgmss(1,ss_name=layer)
>>> Lk.L.Gs.node[1]['ss_name']=layer
>>> Lk.L.g2npy()
>>> Lk.L.save()
>>> fGHz=np.linspace(2,11,181)
>>> #Aa = Antenna('Omni',fGHz=fGHz)
... #Aa = Antenna('Omni',fGHz=fGHz)
... ak,tauk=Lk.eval(force=True)
>>> plt.stem(Lk.H.taud,Lk.H.ak)
>>> plt.stem(Lk.H.taud,Lk.H.ak[:,0,50])
```

```python
>>> Lk.H.ak.shape
(145, 1, 80)
```

```python
>>> cirair = Lk.H.applywavB(wav.sf)

```python
>>> layer = ['PARTITION','PARTITION','PARTITION']
>>> Lk.L.chgmss(1,ss_name=layer)
>>> Lk.L.Gs.node[1]['ss_name']=layer
>>> Lk.L.g2npy()
>>> Lk.L.save()
>>> Lk.eval(force=True)
>>> cirpart = Lk.H.applywavB(wav.sf)
>>> cirpart.plot(typ=['v'],xmin=10,xmax=80)
```

```python
>>> layer = ['METAL','METAL','METAL']
>>> Lk.L.chgmss(1,ss_name=layer)
>>> Lk.L.Gs.node[1]['ss_name']=layer
>>> Lk.L.g2npy()
>>> Lk.L.save()
>>> Lk.eval(force=True)
>>> cirmet = Lk.H.applywavB(wav.sf)
>>> cirmet.plot(typ=['v'],xmin=20,xmax=80)
```

```python
>>> #fig2=plt.figure()
... f,a=cirair.plot(typ=['l20'],color='b')
>>> plt.axis([0,120,-120,-40])
>>> plt.title('A simple illustration of shadowing effect')
>>> plt.legend(['air'])
>>> f,a=cirpart.plot(typ=['l20'],color='k')
>>> plt.axis([0,120,-120,-40])
>>> plt.legend(['wood'])
>>> f,a=cirmet.plot(typ=['l20'],color='r')
>>> plt.axis([0,120,-120,-40])
>>> plt.legend(['metal'])
```

We have modified successively the nature of the 3 surfaces in the sub segment placed in the sepataion partition. The first was AIR, the second WOOD and the third METAL. As the subsegment is placed on the LOS path the blockage effect is clearly visible.
The chosen antennas were omnidirectional `Antenna('Omni')`
