```python
>>> import pylayers.signal.standard as std
```

Parameters of different standards are stored in `wstd.json` file.
The class which handles wireless standards is ` standard`

```python
>>> ws1 = std.Wstandard('ieee80211b')
```

To list the available standard :

```python
>>> ws1.ls()
```

```python
>>> wifi= ws1.load('ieee80211b')
```

```python
>>> ws1
```

```python
>>> AP1 = std.AP()
```

```python
>>> ws2 = std.Wstandard('ieee80211g')
```

```python
>>> ws2
```

```python
>>> ws3 = std.Wstandard('ieee80211a')
```

```python
>>> ws3 = std.Wstandard('bluetooth-class3')
```

```python
>>> ws3
```
