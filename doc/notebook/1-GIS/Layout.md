# Description of the propagation environment

The `Layout` class contains the data structure for describing an Indoor environment. It contathe different graphs helping the implementation of the ray tracing. The class is implemented in the 
[`layout.py`](http://pylayers.github.io/pylayers/modules/pylayers.gis.layout.html)  module.ins

```python
>>> from pylayers.gis.layout import *
>>> from IPython.display import Image
>>> import os
>>> %matplotlib inline
```

## Getting the list of all available Layouts : the `ls()` method

Creating a default Layout is as simple as :

```python
>>> L=Layout()
>>> L

----------------
defstr.ini
----------------

Number of points  : 8
Number of segments  : 9
Number of sub segments  : 3
Number of cycles  : 3
Number of rooms  : 2
degree 0 : []
degree 1 : []
number of node point of degree 2 : 6
number of node point of degree 3 : 2

xrange :(758.49, 768.516)
yrange :(1111.9, 1115.963)

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

Querying the file name associated with the Layout

```python
>>> L.filename
'defstr.ini'
```

 The Layout is described in an `.ini` file.The `ls()` method lists the layout files which are available in the `struc` directory of your current project, which is set up via the $BASENAME environment variable which is crucial to be early defined in order PyLayers find its way to the good directories.

```python
>>> L.ls('ini')
['11Dbibli.ini',
 'CORM1.ini',
 'DLR.ini',
 'DLR2.ini',
 'MADRID-METIS.ini',
 'MOCAP-small.ini',
 'MOCAP-small2.ini',
 'MOCAP-small3.ini',
 'MOCAP.ini',
 'MOCAPext.ini',
 'Scene.ini',
 'TA-Office.ini',
 'TA-OfficeAir.ini',
 'Test_layout6.ini',
 'W2PTIN.ini',
 'WHERE1.ini',
 'WHERE2.ini',
 'd24.ini',
 'defstr.ini',
 'defstr3.ini',
 'homeK_vf.ini',
 'klepal.ini',
 'nicta.ini',
 'rr.ini',
 'scat1.ini',
 'scat2.ini',
 'scattering.ini',
 'test.ini']
```

```python
>>> L=Layout('defstr.ini')
```

```python
>>> L

----------------
defstr.ini
----------------

Number of points  : 8
Number of segments  : 9
Number of sub segments  : 3
Number of cycles  : 3
Number of rooms  : 2
degree 0 : []
degree 1 : []
number of node point of degree 2 : 6
number of node point of degree 3 : 2

xrange :(758.49, 768.516)
yrange :(1111.9, 1115.963)

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
>>> f,a=L.showG('s',nodes=True,slab=True,subseg=True,figsize=(10,10),labels=True)
```

L.ax is  :  (xmin,xmax,ymin,ymax)

```python
>>> L.ax
(758.49, 768.516, 1111.9, 1115.963)
```

```python
>>> L.build()
```

```python
>>> L.ma
(763.49,1115.931)
(768.49,1115.963)
(768.516,1111.964)
(763.516,1111.932)
(758.516,1111.9)
(758.49,1115.9)

vnodes : (-5 6 -4 7 -3 8 -2 9 -1 4 -6 5 )
```

This Layout has 3 cycles. Negative index cycle are outdoor and positive index cycle are indoor. The list of diffraction point for indoor is in ldiffin and the list of diffraction points for outdoor diffraction is in ldiffout. These two listr are

```python
>>> L.Gv.node
{1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 7: {}, 8: {}, 9: {}}
```

```python
>>> L.ldiff
[-6, -4, -3, -1]
```

```python
>>> L.ldiffin
[]
```

```python
>>> L.ldiffout
[-6, -4, -3, -1]
```

```python
>>> L.Gt.node
{-1: {'indoor': False,
  'inter': [(6, -1),
   (6, -1, 0),
   (6, 0, -1),
   (7, -1),
   (7, -1, 0),
   (7, 0, -1),
   (8, -1),
   (8, -1, 0),
   (8, 0, -1),
   (9, -1),
   (9, -1, 0),
   (9, 0, -1),
   (4, -1),
   (4, -1, 0),
   (4, 0, -1),
   (5, -1),
   (5, -1, 0),
   (5, 0, -1),
   (-4,),
   (-3,),
   (-1,),
   (-6,)],
  'polyg': (753.49,1106.9)
  (753.49,1120.963)
  (773.516,1120.963)
  (773.516,1106.9)

  vnodes : (-5 6 -4 7 -3 8 -2 9 -1 4 -6 5 )},
 1: {'cycle': cycle nstr[-8  2 -2  8 -3  7 -4  6 -5  3 -7  1]
  point number 6
  segment number 6
  area : 19.995824
  centroid : [  766.00300113  1113.94747911],
  'indoor': True,
  'inter': [(2, 1),
   (2, 1, 2),
   (2, 2, 1),
   (8, 1),
   (8, 1, 0),
   (8, 0, 1),
   (7, 1),
   (7, 1, 0),
   (7, 0, 1),
   (6, 1),
   (6, 1, 0),
   (6, 0, 1),
   (3, 1),
   (3, 1, 2),
   (3, 2, 1),
   (1, 1),
   (1, 1, 2),
   (1, 2, 1)],
  'isopen': False,
  'polyg': (763.506,1113.432)
  (763.516,1111.932)
  (768.516,1111.964)
  (768.49,1115.963)
  (763.49,1115.931)
  (763.5,1114.432)

  vnodes : (-8 2 -2 8 -3 7 -4 6 -5 3 -7 1 )},
 2: {'cycle': cycle nstr[-8  2 -2  9 -1  4 -6  5 -5  3 -7  1]
  point number 6
  segment number 6
  area : -19.998327
  centroid : [  761.0028967   1113.91576981],
  'indoor': True,
  'inter': [(2, 2),
   (2, 2, 1),
   (2, 1, 2),
   (9, 2),
   (9, 2, 0),
   (9, 0, 2),
   (4, 2),
   (4, 2, 0),
   (4, 0, 2),
   (5, 2),
   (5, 2, 0),
   (5, 0, 2),
   (3, 2),
   (3, 2, 1),
   (3, 1, 2),
   (1, 2),
   (1, 2, 1),
   (1, 1, 2)],
  'isopen': False,
  'polyg': (763.506,1113.432)
  (763.516,1111.932)
  (758.516,1111.9)
  (758.49,1115.9)
  (763.49,1115.931)
  (763.5,1114.432)

  vnodes : (-8 2 -2 9 -1 4 -6 5 -5 3 -7 1 )}}
```

```python
>>> f,a=L.showG
```

```python
>>> L=Layout('DLR.ini')
```

```python
>>> f,a=L.showG('s')
```

To check which are the used slabs :

```python
>>> Slabs = np.unique(L.sla)
>>> for s in Slabs:
...     if s in L.sl:
...         print L.sl[s]
3D_WINDOW_GLASS : GLASS | AIR | GLASS | [0.005, 0.005, 0.005]

AIR : AIR | [0.02]

DOOR : WOOD | [0.03]

METAL : METAL | [0.1]

PARTITION : PLASTER | [0.1]

WALL : BRICK | [0.07]
```

Let's load an other layout

```python
>>> L=Layout('WHERE1.ini')
```

The showG method provides many vizualization of the layout

```python
>>> f,a=L.showG('s',airwalls=False,figsize=(20,10))
```

```python
>>> L=Layout('W2PTIN.ini')
```

```python
>>> f,a = L.showG('s')
```

## The useful numpy arrays of the Layout

The layout data structure is a mix between graph and numpy array. 
numpy arrays are used when high performance is required while graph 
structure is convenient when dealing with different specific tasks. 
The tricky thing for the mind is to have to transcode between node index 
excluding 0 and numpy array index including 0. Below are listed various 
useful numpy array which are mostly used internally.

+ tsg : get segment index in Gs from tahe
+ isss :  sub-segment index above Nsmax
+ tgs : get segment index in tahe from Gs
+ lsss : list of segments with sub-segment
+ sla : list of all slab names (Nsmax+Nss+1)
+ degree : degree of nodes

### `pt` the array of points

The point coordinates are stored in two different places (which in principle  is a bad thing to do !).

    L.Gs.pos : in a dictionnary form (key is the point negative index)
    L.pt : in a numpy array

```python
>>> print np.shape(L.pt)
>>> print len(filter(lambda x: x<0,L.Gs.pos))
```

This dual storage is chosen (temporarily ? ) for computational efficiency reason. The priority goes to the graph and the numpy array is calculated at the end of the edition in the `Layout.g2npy` method (graph to numpy) which is in charge of the conversion.

### tahe (tail-head)

`tahe` is a $(2\times N_{s})$  where $N_s$ denotes the number of segment. The first line  is the tail index of the segment $k$ and the second line is the head of the segment $k$. Where $k$ is the index of a given segment (starting in 0).

```python
>>> L.build()
```

The figure below illustrates a Layout and a surimposition of the graph of cycles $\mathcal{G}_c$. Those cycles are automatically extracted from a well defined layout. This concept of **cycles** is central in the ray determination algorithm which is implemented in PyLayers. Notice that the exterior region is the cycle indexed by 0. All the rooms which have a common frontier with the exterior cycle are here connected to the origin (corresponding to exterior cycle).

```python
>>> f,a = L.showG('s')
>>> nx.draw(L.Gc,L.Gc.pos)
```

```python
>>> nx.draw_networkx_nodes(L.Gi,L.Gi.pos,node_color='blue',node_size=1)
>>> nx.draw_networkx_edges(L.Gi,L.Gi.pos,node_color='blue',node_size=1)
```

## `tgs` : trancodage from graph indexing to numpy array indexing

`tgs` is an array with length $N_s$+1. The index 0 is not used because none segment has 0 as an index.

```python
>>> ns = 5
>>> utahe = L.tgs[ns]
```

```python
>>> tahe =  L.tahe[:,utahe]
```

```python
>>> ptail = L.pt[:,tahe[0]]
>>> phead = L.pt[:,tahe[1]]
```

```python
>>> print ptail
```

```python
>>> print phead
```

```python
>>> L.Gs.node[5]
```

```python
>>> print L.Gs.pos[-8]
>>> print L.Gs.pos[-139]
```

```python
>>> aseg = np.array([4,7,134])
```

```python
>>> print np.shape(aseg)
```

```python
>>> pt  = L.tahe[:,L.tgs[aseg]][0,:]
>>> ph = L.tahe[:,L.tgs[aseg]][1,:]
>>> pth = np.vstack((pt,ph))
```

```python
>>> np.shape(pth)
```

## `Layout.seg2pts` a function for getting points coordinates from segment number array

```python
>>> L.seg2pts(aseg)
```

```python
>>> aseg = array(filter(lambda x: x>0,L.Gs.nodes()))
>>> pth = L.seg2pts(aseg)
```

```python
>>> from pylayers.util.plotutil import displot
```

```python
>>> displot(pth[0:2,:],pth[2:,:])
>>> plt.axis('off')
```
