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
```

Let start by loading a simple layout with 2 single rooms. The multi subsegment appears in the middle with the red vertical lines. Each subsegment is materialized by a  segment.

```python
>>> L=Layout('defstr3.ini')
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
structure saved in  defstr3.str2
structure saved in  defstr3.ini
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
>>> tx=np.array([759,1114,1.5])
>>> rx=np.array([767,1114,1.5])
```

```python
>>> L.chgmss(1,ss_name=['WOOD','AIR','WOOD'],ss_z =[(0.0,2.7),(2.7,2.8),(2.8,3)],ss_offset=[0,0,0])
>>> L.save()
>>> fGHz=np.linspace(1,11,100)
>>> #Aa = Antenna('S1R1.vsh3')
... #Ab = Antenna('S1R1.vsh3')
... Aa = Antenna('Omni',fGHz=fGHz)
>>> Ab = Antenna('Omni',fGHz=fGHz)
>>> Lk = DLink(L=L,a=tx,b=rx,Aa=Aa,Ab=Ab,fGHz=np.linspace(1,11,100))
structure saved in  defstr3.str2
structure saved in  defstr3.ini
```

A link is the set of a layout and 2 termination points.

```python
>>> f,a=Lk.show()
```

On the figure above, we can see the Tx and Rx each placed in a different room appart from a wall with a subsegement placed in the middle.
Then for evaluating the radio link, simply type:

```python
>>> ak,tauk=Lk.eval(force=True,a=tx,b=rx)
checkh5
Start Signatures
algo 7
Signatures'> from 2_1_3 saved
Stop signature 0.0677840709686
Start Rays
Rays'> from 3_2_1 saved
Stop rays 0.322803974152
Ctilde'> from 2_1_0 saved
```

At that point the channel has been evaluated and all the data stored in an `hdf5` file

## Link members

The Signature of the radio channel is in `Lk.Si`, the 3D rays are in `Lk.R`, the propagation channel is in `Lk.C` and the transmission channel is in `Lk.H`

```python
>>> Lk.R
Rays3D
----------
1 / 1 : [0]
2 / 6 : [1 2 3 4 5 6]
3 / 19 : [ 7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]
4 / 40 : [26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65]
5 / 49 : [ 66  67  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  83
  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 101
 102 103 104 105 106 107 108 109 110 111 112 113 114]
6 / 34 : [115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132
 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148]
7 / 6 : [149 150 151 152 153 154]
-----
ni : 721
nl : 1597
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
structure saved in  defstr3.str2
structure saved in  defstr3.ini
checkh5
Start Signatures
algo 7
Stop signature 0.0449650287628
Signatures'> from 2_1_3 saved
Start Rays
Rays'> from 3_2_1 saved
Stop rays 0.557983160019
Ctilde'> from 2_1_0 saved
Tchannel'> from 2_1_0_0_0_2_2 saved
<Container object of 3 artists>
```

```python
>>> cirair = Lk.H.applywavB(wav.sf)
```

```python
>>> cirair.plot(typ=['v'],xmin=20,xmax=80)
(<matplotlib.figure.Figure at 0x7fad347f2ad0>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fad347f9ad0>]], dtype=object))
```

```python
>>> layer = ['PARTITION','PARTITION','PARTITION']
>>> Lk.L.chgmss(1,ss_name=layer)
>>> Lk.L.Gs.node[1]['ss_name']=layer
>>> Lk.L.g2npy()
>>> Lk.L.save()
>>> Lk.eval(force=True)
>>> cirpart = Lk.H.applywavB(wav.sf)
>>> cirpart.plot(typ=['v'],xmin=20,xmax=80)
structure saved in  defstr3.str2
structure saved in  defstr3.ini
checkh5
Start Signatures
algo 7
Stop signature 0.0465669631958
Signatures'> from 2_1_3 saved
Start Rays
Rays'> from 3_2_1 saved
Stop rays 0.558837890625
Ctilde'> from 2_1_0 saved
Tchannel'> from 2_1_0_0_0_2_2 saved
(<matplotlib.figure.Figure at 0x7fad2f51ab10>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fad3463c7d0>]], dtype=object))
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
structure saved in  defstr3.str2
structure saved in  defstr3.ini
checkh5
Start Signatures
algo 7
Stop signature 0.0461058616638
Signatures'> from 2_1_3 saved
Start Rays
Rays'> from 3_2_1 saved
Stop rays 0.56045293808
Ctilde'> from 2_1_0 saved
Tchannel'> from 2_1_0_0_0_2_2 saved
(<matplotlib.figure.Figure at 0x7fad3629f190>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fad3465b790>]], dtype=object))
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
