Robust Geometric Positioning Algorithm

```python
>>> from pylayers.location.geometric.constraints.cla import *
>>> from pylayers.location.geometric.constraints.toa import *
>>> from pylayers.location.geometric.constraints.tdoa import *
>>> from pylayers.location.geometric.constraints.rss import *
>>> from pylayers.network.model import *
>>> 
>>> import matplotlib.pyplot as plt
```

Let's define 4 anchors in the plane.

```python
>>> pt1=np.array(([0,0]))
>>> pt2=np.array(([10,15]))
>>> pt3=np.array(([5,28]))
>>> pt4=np.array(([-10,-10]))
>>> # the unknown point is p
... p = np.array((0,5))
```

Then displaying the scene with ancho nodes (in red) and blind nodes in blue

```python
>>> f=plt.figure()
>>> ax=f.add_subplot(111)
>>> ax.plot(pt1[0],pt1[1],'or',label='anchor 1')
>>> ax.plot(pt2[0],pt2[1],'or',label='anchor 2')
>>> ax.plot(pt3[0],pt3[1],'or',label='anchor 3')
>>> ax.plot(pt4[0],pt4[1],'or',label='anchor 4')
>>> ax.plot(p[0],p[1],'xb',label='p')
>>> 
>>> 
>>> ax.legend(loc='best')
```

The euclidian distance betwenn blind node $p$ and anchors ans the corresponding Time of flight are evaluated.

```python
>>> d1=np.sqrt(np.sum((pt1-p)**2))
>>> d2=np.sqrt(np.sum((pt2-p)**2))
>>> d3=np.sqrt(np.sum((pt3-p)**2))
>>> d4=np.sqrt(np.sum((pt4-p)**2))
>>> 
>>> toa1=d1/0.3
>>> toa2=d2/0.3
>>> toa3=d3/0.3
>>> toa4=d4/0.3
>>> 
>>> print 'distance p-1=',d1, '/ toa1=',toa1
>>> print 'distance p-2=',d2, '/ toa2=',toa2
>>> print 'distance p-3=',d3, '/ toa3=',toa3
>>> print 'distance p-4=',d4, '/ toa3=',toa4
distance p-1= 5.0 / toa1= 16.6666666667
distance p-2= 14.1421356237 / toa2= 47.1404520791
distance p-3= 23.5372045919 / toa3= 78.4573486396
distance p-4= 18.0277563773 / toa3= 60.0925212577
```

# RGPA (Robust Geometric Positioning Algorithm)

## Exploiting Time of Arrival (ToA) only

We call `Constraint Layer Array` (CLA) the object which gathers all the geometric constraints of a considered scenario.

```python
>>> C=CLA()
```

Instanciate TOA constraints, notice that their id are differents

```python
>>> T1=TOA(id=0,value = toa1, std = np.array([1.0]), p = pt1)
>>> T2=TOA(id=1,value = toa2, std = np.array([1.0]), p = pt2)
>>> T3=TOA(id=2,value = toa3, std = np.array([1.0]), p = pt3)
>>> T4=TOA(id=3,value = toa4, std = np.array([1.0]), p = pt4)
```

Add TOA contraints to the CLA

```python
>>> C.append(T1)
>>> C.append(T2)
>>> C.append(T3)
>>> C.append(T4)
```

All the constraints of the CLA can be listed as follows

```python
>>> C.c
[node | peer   |type | rat  | p              | value    | std  | runable| usable|
    0 |        |TOA  |      | [0 0]          | [ 16.667]| [ 1.]|       1|      1|,
 node | peer   |type | rat  | p              | value    | std  | runable| usable|
    1 |        |TOA  |      | [10 15]        | [ 47.14] | [ 1.]|       1|      1|,
 node | peer   |type | rat  | p              | value    | std  | runable| usable|
    2 |        |TOA  |      | [ 5 28]        | [ 78.457]| [ 1.]|       1|      1|,
 node | peer   |type | rat  | p              | value    | std  | runable| usable|
    3 |        |TOA  |      | [-10 -10]      | [ 60.093]| [ 1.]|       1|      1|]
```

Get information on the cla :

  + type : TOA / RSS
  + p : Position of the origin of the constraint
  + value : power ( RSS ) / time in ns ( TOA)
  + std : standard deviation $\sigma^2$ of value
  + runable : does the constraint has a position p ?
  + obsolete : does the value has been obtained recently ?
  + usuable : runbale AND NOT obsolete
  + evlauated : obsolete

```python
>>> C.info()
type , p              , value, std  , runable, usable, obsolete, evaluated
TOA  , [0 0]          , [ 16.667], [ 1.],       1,      1,        0,         0
type , p              , value, std  , runable, usable, obsolete, evaluated
TOA  , [10 15]        , [ 47.14], [ 1.],       1,      1,        0,         0
type , p              , value, std  , runable, usable, obsolete, evaluated
TOA  , [ 5 28]        , [ 78.457], [ 1.],       1,      1,        0,         0
type , p              , value, std  , runable, usable, obsolete, evaluated
TOA  , [-10 -10]      , [ 60.093], [ 1.],       1,      1,        0,         0
```

Update the CLA

```python
>>> C.update()
```

Compute the cla

```python
>>> C.compute()
cluster
enter in HT processing
!!!!! HT FAIL !!!!!!!
2 first constraint of CLA have to be TOA and others RSS in order to use HT
/home/uguen/Documents/rch/devel/pylayers/pylayers/location/geometric/util/boxn.py:92: FutureWarning: comparison to `None` will result in an elementwise object comparison in the future.
  if Lb==None:

True
```

show the estimated position

```python
>>> C.pe
array([ -4.735e-03,   4.992e+00])
```

to be compare with the actual position value

```python
>>> p
array([0, 5])
```

## RSS

The RSS is a quantity which is weakly related to distance via a parametric model. The bettet the model, better would be tthe inference ab out tthe associated distance.
t
To model the Path Loss shadowing model is widely used.

To define the classical path loss shadowing model widely used in this context the `PLSmodel` class has been defined.

```python
>>> M = PLSmodel(f=3.0,rssnp=2.64,d0=1.0,sigrss=3.0,method='mode')
```

For simulation purpose : get RSS from distances (or associated delay) with the above model

```python
>>> toa1
16.666666666666668
```

```python
>>> M.getPL(toa1,1)
10.393498107059024
```

## TDOA

```python
>>> Td1=TDOA(id=0,value = toa1-toa2, std = np.array([1.0]), p = np.array([pt1,pt2]))
>>> Td2=TDOA(id=1,value = toa1-toa3, std = np.array([1.0]), p = np.array([pt1,pt3]))
>>> Td3=TDOA(id=2,value = toa1-toa4, std = np.array([1.0]), p = np.array([pt1,pt4]))
```

```python
>>> C=CLA()
>>> C.append(Td1)
>>> C.append(Td2)
>>> C.append(Td3)
```

```python
>>> C.compute()
TDOA 2.0
TDOA 2.0
TDOA 2.0
TDOA 1.5
TDOA 1.5
TDOA 1.5
TDOA 1.375
TDOA 1.375
TDOA 1.375
TDOA 1.375
TDOA 1.375
TDOA 1.375
True
```

```python
>>> C.pe
array([ 0.021,  4.987])
```

```python
>>> from IPython.core.display import HTML
>>> 
>>> def css_styling():
...     styles = open("../styles/custom.css", "r").read()
...     return HTML(styles)
>>> css_styling()
<IPython.core.display.HTML at 0x7fc7801770d0>
```
