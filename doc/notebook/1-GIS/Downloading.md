```python
>>> from pylayers.gis.ezone import *
>>> %matplotlib inline
```

An Ezone is an objet which gathers vector and raster information from a 1 degree by 1 degree tile of earth

```python
>>> prefix=enctile(-1.5,10.5)
>>> print prefix
N10W002
```

```python
>>> prefix
'N10W002'
```

```python
>>> dectile(prefix=prefix)
(-2, -1, 10, 11)
```

```python
>>> int('08')
8
```

```python
>>> E=Ezone(prefix)
```

```python
>>> E.prefix
'N10W002'
```

An `Ezone` can be obtained from a point (longitude,Latitude)

```python
>>> r=E.getdem()
```

```python
>>> E.saveh5()
```

```python
>>> f,a = E.show(source='srtm',clim=[0,500])
```

```python
>>> from IPython.core.display import HTML
>>> 
```

