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
Download srtm file
SRTMDownloader - server= dds.cr.usgs.gov, directory=srtm/version2_1/SRTM3.
no aster file for this point
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
>>> def css_styling():
...     styles = open("../styles/custom.css", "r").read()
...     return HTML(styles)
>>> css_styling()
<IPython.core.display.HTML at 0x7f7207397690>
```
