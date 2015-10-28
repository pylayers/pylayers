```python
>>> import numpy as np 
>>> import scipy as sp
>>> import networkx as nx
>>> %pylab inline
Populating the interactive namespace from numpy and matplotlib
```

```python
>>> pt = 10*np.random.rand(3,50)
>>> tx = 10*np.random.rand(3,1)
>>> rx = 10*np.random.rand(3,1)
>>> CardEt = 6
>>> CardEr = 6
>>> CardEs = 40
```

```python
>>> taue = np.sqrt(np.sum((tx-rx)*(tx-rx),axis=0))/0.3
```

```python
>>> dpttx = np.sqrt(np.sum((pt-tx)*(pt-tx),axis=0))
>>> u = np.argsort(dpttx)
>>> dptrx = np.sqrt(np.sum((pt-rx)*(pt-rx),axis=0))
>>> v = np.argsort(dptrx)
```

```python
>>> u
array([ 5, 11, 10, 20, 22,  4, 29, 24, 17, 35, 31, 19, 47, 27, 28, 41,  2,
       26, 44, 14, 33, 48, 46, 36, 40,  3,  6, 18, 21, 39, 42, 37,  8, 15,
        9, 25,  1, 30,  0, 49, 13, 43, 45, 12, 16, 34,  7, 38, 32, 23])
```

```python
>>> v
array([ 9, 37, 34,  1,  8, 36, 48, 39, 26, 49, 14,  3, 41, 16, 47, 15, 35,
       19, 31, 30,  2, 10, 29, 45, 18, 27,  0, 38, 21, 32, 20, 42, 24,  4,
       40, 28, 22,  5, 44,  7, 11, 43, 23, 25, 12, 17, 13, 33, 46,  6])
```

```python
>>> fGHz=2.4
```

```python
>>> Ed = 1/(4*np.pi*fGHz*taue)**2
```

```python
>>> Ed
array([  1.01027240e-06])
```

```python
>>> G = nx.Graph()
>>> G.pos ={}
>>> G.add_node(-1)
>>> G.pos[-1]=(tx[0],tx[1])
>>> G.add_node(-2)
>>> G.pos[-2]=(rx[0],rx[1])
```

```python
>>> for k in range(CardEt):
...     G.add_node(u[k])
...     G.pos[u[k]]=(pt[0,u[k]],pt[1,u[k]])
...     G.add_edge(-1,u[k])
>>> 
>>> for k in range(CardEr):
...     G.add_node(v[k])
...     G.pos[v[k]]=(pt[0,v[k]],pt[1,v[k]])
...     G.add_edge(-2,v[k])
...     
>>> for k in range(CardEs):
...     n1 = np.random.randint(1,50,1)[0]
...     n2 = np.random.randint(1,50,1)[0]
...     if n1!=n2:
...         if n1 not in G.node:
...             G.add_node(n1)
...             G.pos[n1]=(pt[0,n1],pt[1,n1])
...         if n2 not in G.node:
...             G.add_node(n2)
...             G.pos[n2]=(pt[0,n2],pt[1,n1])
...         G.add_edge(n1,n2)
```

```python
>>> G
<networkx.classes.graph.Graph at 0x7f7251c0c210>
```

```python
>>> nx.draw(G,G.pos)
```
