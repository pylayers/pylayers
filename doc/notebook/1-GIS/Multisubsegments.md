# Effect of Modyfiying the Nature of Sub-Segments

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

This section presents a simple Ray Tracing simulation with different material properties of a subsegment separating 2 rooms.

```python
>>> L=Layout('defstr3.ini')
>>> f,a=L.showG('s',subseg=True,figsize=(10,10))
```

The studied configuration is composed of a simple 2 rooms building separated by a subsegment which has a multi subsegment attribute. The attribute of the subsegment can be changed  with the method [`chgmss`](http://pylayers.github.io/pylayers/modules/generated/pylayers.gis.layout.Layout.chgmss.html) (change multisubsegment)

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
 9: {'connect': [-1, -2],
  'name': 'WALL',
  'ncycles': [2, 0],
  'norm': array([-0.00639987,  0.99997952,  0.        ]),
  'offset': 0,
  'transition': False,
  'z': (0.0, 3.0)}}
```

We define now two points which are the termination of a radio link.

```python
>>> tx=np.array([759,1114,1.0])
>>> rx=np.array([767,1114,1.5])
```

```python
>>> L.chgmss(1,ss_name=['WOOD','AIR','WOOD'],ss_z =[(0.0,2.7),(2.7,2.8),(2.8,3)],ss_offset=[0,0,0])
>>> L.save()
>>> fGHz=np.linspace(1,11,100)
>>> Aa = Antenna('S1R1.vsh3',fGHz=fGHz)
>>> Ab = Antenna('S1R1.vsh3',fGHz=fGHz)
>>> Lk = DLink(L=L,a=tx,b=rx,Aa=Aa,Ab=Ab,fGHz=np.linspace(1,11,100))
structure saved in  defstr3.str2
structure saved in  defstr3.ini
```

```python
>>> f,a=Lk.show()
```

On the figure above, we can see the Tx and Rx each placed in a different room appart from a wall with a subsegement placed in the middle.
Then for evaluating the radio link, simply type:

```python
>>> ak,tauk=Lk.eval(force=True)
Signatures'> from 2_1_3 saved
Rays'> from 3_0_1 saved
Ctilde'> from 0_1_0 saved
Tchannel'> from 0_1_0_0_0_1_1 saved
```

At that point the channel has been evaluated and all the data stored in an `hdf5` file

## The different members of the link are

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
>>> f,a=Lk.C.show(cmap='jet',fig=f,typ='l20')
```

```python
>>> Lk.H
freq : 1.0 5.95 100
shape  : (155, 100)
tau (min, max) : 26.7186992365 95.816242908
dist :8.01560977094 28.7448728724
Friis factor -j c/(4 pi f) has been applied
 calibrated : No
 windowed : No
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
>>> Lk.H.taud
array([ 26.71869924,  27.93842436,  29.10708199,  29.64889324,
        30.03048589,  30.075433  ,  36.72255959,  30.75261837,
        31.12068041,  31.1640552 ,  31.81807982,  32.17395555,
        32.21591227,  32.36081306,  32.66533294,  33.05244038,
        34.35921355,  37.6193175 ,  37.62655033,  37.86521462,
        38.49519081,  38.90590625,  39.23286412,  40.03698677,
        60.02330145,  60.0237882 ,  33.67032288,  34.04600538,
        34.64617309,  34.81939407,  35.01138599,  35.14489175,
        35.18330575,  36.68410469,  36.99319869,  37.02969534,
        38.50225911,  38.7355282 ,  39.35848858,  39.58671257,
        39.75344272,  40.04709276,  40.07348617,  40.27031144,
        40.58328319,  40.8610692 ,  40.89683313,  41.01072684,
        41.66885699,  42.04868734,  42.35118067,  42.60531715,
        45.28154414,  45.47942049,  47.70172314,  48.49380416,
        60.57609581,  60.57657812,  61.12389099,  61.12436898,
        61.38177086,  61.38655397,  61.589339  ,  61.5932084 ,
        63.33967443,  70.03865417,  37.42134832,  37.75972919,
        39.16236258,  39.48582634,  40.87097143,  41.08971465,
        41.6785673 ,  41.82213078,  41.8930939 ,  42.03698148,
        42.83408426,  42.97677133,  43.13106967,  43.27298187,
        43.38690997,  43.59404942,  43.6053373 ,  43.89710511,
        44.00333673,  44.50096862,  44.78709963,  45.49315307,
        46.01179095,  46.20654017,  46.73062778,  46.9223936 ,
        48.00940863,  48.19516497,  48.39546526,  49.07940224,
        49.17637348,  49.84959754,  61.92243907,  61.92718042,
        62.12820088,  62.13203673,  62.45842719,  62.46312786,
        62.66242903,  62.66623217,  62.73858502,  62.73905071,
        63.7923458 ,  63.7928038 ,  63.86376926,  64.38359799,
        70.51297572,  70.98412788,  93.34926748,  44.01253199,
        44.21573608,  45.50204727,  45.69862854,  45.84130714,
        46.11893144,  47.27323528,  47.54249857,  48.69876778,
        48.82367841,  48.88190455,  49.00725478,  49.37850393,
        49.55912892,  50.17055816,  50.34922397,  51.07629317,
        51.81681556,  52.36526576,  53.0878113 ,  64.039481  ,
        64.04406563,  64.23846209,  64.24217193,  65.07217885,
        65.07669072,  65.26801165,  65.27166298,  65.91849278,
        66.92220128,  72.37918493,  73.29447281,  93.70566901,
        94.06072013,  51.36376787,  51.53743551,  52.64570242,
        52.81515495,  95.11792193,  95.81624291])
```

```python
>>> cir = Lk.H.applywavB(wav.sf)
> /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/channel.py(2510)ft1()
-> if len(tau) == 1:
(Pdb) p r
TUsignal :  (1239,)  (1, 1239) 
(Pdb) p x
array([-9.99192897, -9.97578692, -9.95964487, ...,  9.95964487,
        9.97578692,  9.99192897])
(Pdb) n
> /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/channel.py(2513)ft1()
-> for i in range(len(tau)):
(Pdb) p len(tau)
155
(Pdb) n
> /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/channel.py(2514)ft1()
-> si = bs.TUsignal(self.s.x, self.s.y[i, :])
(Pdb) n
> /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/channel.py(2515)ft1()
-> si.translate(tau[i])
(Pdb) p si.x
array([-9.99192897, -9.97578692, -9.95964487, ...,  9.95964487,
        9.97578692,  9.99192897])
(Pdb) p tau[i]
26.718699236469
(Pdb) n
> /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/channel.py(2516)ft1()
-> r = r + si
(Pdb) p si.x
array([ 16.72677026,  16.74291231,  16.75905436, ...,  36.67834411,
        36.69448616,  36.71062821])
(Pdb) n
> /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/channel.py(2513)ft1()
-> for i in range(len(tau)):
(Pdb) p r
TUsignal :  (1239,)  (1, 1239) 
(Pdb) p r.x
array([-9.99192897, -9.97578692, -9.95964487, ...,  9.95964487,
        9.97578692,  9.99192897])
```

```python
>>> Lk.H.taud
array([ 26.71869924,  27.93842436,  29.10708199,  29.64889324,
        30.03048589,  30.075433  ,  36.72255959,  30.75261837,
        31.12068041,  31.1640552 ,  31.81807982,  32.17395555,
        32.21591227,  32.36081306,  32.66533294,  33.05244038,
        34.35921355,  37.6193175 ,  37.62655033,  37.86521462,
        38.49519081,  38.90590625,  39.23286412,  40.03698677,
        60.02330145,  60.0237882 ,  33.67032288,  34.04600538,
        34.64617309,  34.81939407,  35.01138599,  35.14489175,
        35.18330575,  36.68410469,  36.99319869,  37.02969534,
        38.50225911,  38.7355282 ,  39.35848858,  39.58671257,
        39.75344272,  40.04709276,  40.07348617,  40.27031144,
        40.58328319,  40.8610692 ,  40.89683313,  41.01072684,
        41.66885699,  42.04868734,  42.35118067,  42.60531715,
        45.28154414,  45.47942049,  47.70172314,  48.49380416,
        60.57609581,  60.57657812,  61.12389099,  61.12436898,
        61.38177086,  61.38655397,  61.589339  ,  61.5932084 ,
        63.33967443,  70.03865417,  37.42134832,  37.75972919,
        39.16236258,  39.48582634,  40.87097143,  41.08971465,
        41.6785673 ,  41.82213078,  41.8930939 ,  42.03698148,
        42.83408426,  42.97677133,  43.13106967,  43.27298187,
        43.38690997,  43.59404942,  43.6053373 ,  43.89710511,
        44.00333673,  44.50096862,  44.78709963,  45.49315307,
        46.01179095,  46.20654017,  46.73062778,  46.9223936 ,
        48.00940863,  48.19516497,  48.39546526,  49.07940224,
        49.17637348,  49.84959754,  61.92243907,  61.92718042,
        62.12820088,  62.13203673,  62.45842719,  62.46312786,
        62.66242903,  62.66623217,  62.73858502,  62.73905071,
        63.7923458 ,  63.7928038 ,  63.86376926,  64.38359799,
        70.51297572,  70.98412788,  93.34926748,  44.01253199,
        44.21573608,  45.50204727,  45.69862854,  45.84130714,
        46.11893144,  47.27323528,  47.54249857,  48.69876778,
        48.82367841,  48.88190455,  49.00725478,  49.37850393,
        49.55912892,  50.17055816,  50.34922397,  51.07629317,
        51.81681556,  52.36526576,  53.0878113 ,  64.039481  ,
        64.04406563,  64.23846209,  64.24217193,  65.07217885,
        65.07669072,  65.26801165,  65.27166298,  65.91849278,
        66.92220128,  72.37918493,  73.29447281,  93.70566901,
        94.06072013,  51.36376787,  51.53743551,  52.64570242,
        52.81515495,  95.11792193,  95.81624291])
```

```python
>>> cir.x
array([-9.99192897, -9.97578692, -9.95964487, ...,  9.95964487,
        9.97578692,  9.99192897])
```

```python
>>> layer = ['AIR','AIR','AIR']
>>> L.chgmss(1,ss_name=layer)
>>> L.Gs.node[1]['ss_name']=layer
>>> L.g2npy()
>>> L.save()
>>> Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))
>>> Lk.eval(force=True)
>>> cirair = Lk.H.applywavB(wav.sf)
>>> cirair.plot(typ=['v'],xmin=20,xmax=80)
structure saved in  defstr3.str2
structure saved in  defstr3.ini
Signatures'> from 2_1_3 saved
Rays'> from 3_0_1 saved
Ctilde'> from 0_1_0 saved
Tchannel'> from 0_1_0_0_0_1_1 saved
(<matplotlib.figure.Figure at 0x7f39272475d0>,
 array([[<matplotlib.axes.AxesSubplot object at 0x7f3926fbb690>]], dtype=object))
```

```python
>>> layer = ['PARTITION','PARTITION','PARTITION']
>>> L.chgmss(1,ss_name=layer)
>>> L.Gs.node[1]['ss_name']=layer
>>> L.g2npy()
>>> L.save()
>>> Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))
>>> Lk.eval(force=True)
>>> cirpart = Lk.H.applywavB(wav.sf)
>>> cirpart.plot(typ=['v'],xmin=20,xmax=80)
structure saved in  defstr3.str2
structure saved in  defstr3.ini
Signatures'> from 2_1_3 saved
Rays'> from 3_0_1 saved
Ctilde'> from 0_1_0 saved
Tchannel'> from 0_1_0_0_0_1_1 saved
(<matplotlib.figure.Figure at 0x7f392b898d90>,
 array([[<matplotlib.axes.AxesSubplot object at 0x7f392b8988d0>]], dtype=object))
```

```python
>>> layer = ['METAL','METAL','METAL']
>>> L.chgmss(1,ss_name=layer)
>>> L.Gs.node[1]['ss_name']=layer
>>> L.g2npy()
>>> L.save()
>>> Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))
>>> Lk.eval(force=True)
>>> cirmet = Lk.H.applywavB(wav.sf)
>>> cirmet.plot(typ=['v'],xmin=20,xmax=80)
structure saved in  defstr3.str2
structure saved in  defstr3.ini
Signatures'> from 2_1_3 saved
Rays'> from 3_0_1 saved
Ctilde'> from 0_1_0 saved
Tchannel'> from 0_1_0_0_0_1_1 saved
(<matplotlib.figure.Figure at 0x7f392793d690>,
 array([[<matplotlib.axes.AxesSubplot object at 0x7f3926fb4d10>]], dtype=object))
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
