# Algebraic Localization Class

For localization in the plane it is important not to provide the z coordinates. Otherwise
a singularity problem arises.

```python
>>> import numpy as np
>>> from pylayers.location.algebraic.algebraic  import *
```

```python
>>> c = 0.2997924583
>>> nodes={}
>>> ldp={}
>>> p1 = np.array([-1,-1])
>>> p2 = np.array([1,-1])
>>> p3 = np.array([0,1])
>>> nodes['BN'] = np.array([[0],[0]])
>>> nodes['RN_TOA']=np.vstack((p1,p2,p3)).T
>>> ldp['TOA']=np.array([np.sqrt(2),np.sqrt(2),1])
>>> ldp['TOA_std']=np.array([1,1,1])
>>> Alg = algloc(nodes,ldp)
```

```python
>>> print ldp
>>> Alg.plot()
{'TOA_std': array([1, 1, 1]), 'TOA': array([ 1.41421356,  1.41421356,  1.        ])}
(<matplotlib.figure.Figure at 0x7f9a02ad8a90>,
 <mpl_toolkits.mplot3d.axes3d.Axes3D at 0x7f9a02ad8e50>)
```

```python
>>> pest = Alg.ls_locate(toa=True,tdoa=False,rss=False)
```

```python
>>> Alg.plot()
(<matplotlib.figure.Figure at 0x7f9a02b1b790>,
 <mpl_toolkits.mplot3d.axes3d.Axes3D at 0x7f9a02399a90>)
```

```python
>>> pest
array([[ 0.        ],
       [-0.22753112]])
```

```python
>>> Alg
Nodes : {'BN': array([[0],
       [0]]), 'RN_TOA': array([[-1,  1,  0],
       [-1, -1,  1]])}
LDPs :{'TOA_std': array([1, 1, 1]), 'TOA': array([ 1.41421356,  1.41421356,  1.        ])}
```
