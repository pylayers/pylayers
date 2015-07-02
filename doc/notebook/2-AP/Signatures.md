# Ray Signatures

Signature are calculated from a cycle to an other cycle of the Layout. They are used for the determination of rays once the transmitter and receiver are known. The best algorithmic manner to get a signature is not stabilized yet, a lot of approaches have been implemented so far, with very different performances. It is expected that the signature is delivering all its utility when dealing with mobile trajectories.

The evaluation of a signature from one cycle to another is implemented in the `pylayers.simul.Link.DLink` class.

```python
>>> import time
>>> from pylayers.gis.layout import *
>>> from pylayers.antprop.signature import *
>>> from pylayers.antprop.rays import *
>>> %matplotlib inline
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
```

```python
>>> L = Layout('defstr3.ini')
>>> L.build()
>>> L.dumpw()
```

Showing the graph of rooms with 2 rooms separated by a DOOR segment

```python
>>> L.showG('sv',figsize=(8,8))
>>> a=plt.axis('off')
```

The graph of interactions is shown below.

```python
>>> L.showG('si',figsize=(10,5))
>>> a=plt.axis('off')
```

All the interactions of a given cycle are stored as meta information in nodes of `Gt`

```python
>>> L.Gi.node
{(1, 1): {},
 (1, 1, 2): {},
 (1, 2): {},
 (1, 2, 1): {},
 (2, 1): {},
 (2, 1, 2): {},
 (2, 2): {},
 (2, 2, 1): {},
 (3, 1): {},
 (3, 1, 2): {},
 (3, 2): {},
 (3, 2, 1): {},
 (4, 0): {},
 (4, 0, 2): {},
 (4, 2): {},
 (4, 2, 0): {},
 (5, 0): {},
 (5, 0, 2): {},
 (5, 2): {},
 (5, 2, 0): {},
 (6, 0): {},
 (6, 0, 1): {},
 (6, 1): {},
 (6, 1, 0): {},
 (7, 0): {},
 (7, 0, 1): {},
 (7, 1): {},
 (7, 1, 0): {},
 (8, 0): {},
 (8, 0, 1): {},
 (8, 1): {},
 (8, 1, 0): {},
 (9, 0): {},
 (9, 0, 2): {},
 (9, 2): {},
 (9, 2, 0): {}}
```

```python
>>> L.Gt.node[0]['inter']
[(6, 0),
 (6, 0, 1),
 (6, 1, 0),
 (7, 0),
 (7, 0, 1),
 (7, 1, 0),
 (8, 0),
 (8, 0, 1),
 (8, 1, 0),
 (9, 0),
 (9, 0, 2),
 (9, 2, 0),
 (4, 0),
 (4, 0, 2),
 (4, 2, 0),
 (5, 0),
 (5, 0, 2),
 (5, 2, 0),
 (-4,),
 (-3,),
 (-1,),
 (-6,)]
```

The signature is calculated with as parameters the  Layout object and two cycle numbers.
In example below it is 0 and 1.

```python
>>> Si = Signatures(L,0,1)
```

The cold start determination of the signature is done with a `run` function. The code is not in its final shape here and there is room for significant acceleration in incorporating propagation based heuristics. The mitigation of graph exploration depth is done in setting a `cutoff` value which limits the exploration in the interaction graph.

```python
>>> Si.run5(cutoff=5,diffraction=False,algo='old')
```

The representation method of a signature gives informations about the different signatures. Signatures are grouped by number of interactions.

```python
>>> L.Gt.pos
{0: (758.49, 1111.9),
 1: (766.00300113353387, 1113.947479109665),
 2: (761.00289669547806, 1113.915769812613)}
```

```python
>>> ptx = np.array(L.Gt.pos[0])+np.random.rand(2)
>>> prx = np.array(L.Gt.pos[1])+np.random.rand(2)
>>> print ptx
>>> print prx
[  759.3213136   1112.39386191]
[  766.8885467   1113.98892221]
```
