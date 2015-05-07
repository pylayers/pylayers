# Troels Perdersen Propagation graph

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
array([ 6, 41, 49,  8, 37, 32, 30, 38,  2,  9, 20, 11, 27, 14, 31, 26, 24,
       48, 28, 18,  5, 46, 35, 21,  7, 40,  0, 19,  1, 23, 25, 12, 13,  3,
       45, 17, 47, 34,  4, 15, 42, 44, 36, 29, 16, 33, 39, 22, 10, 43])
```

```python
>>> v
array([19, 24, 15, 26, 36, 30, 12, 21,  8, 32, 42, 38, 28, 14, 45, 49,  0,
       27, 34,  4, 13, 31, 23, 37, 40, 46, 33, 20,  6,  3, 22, 17,  2,  7,
       47, 48,  1, 41, 25,  5, 10, 16,  9, 18, 11, 29, 39, 35, 43, 44])
```

```python
>>> fGHz=2.4
```

```python
>>> Ed = 1/(4*np.pi*fGHz*taue)**2
```

```python
>>> Ed
array([  2.24323622e-06])
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
<networkx.classes.graph.Graph at 0x7f7251a57290>
```

```python
>>> nx.draw(G,G.pos)
```

```python

```
