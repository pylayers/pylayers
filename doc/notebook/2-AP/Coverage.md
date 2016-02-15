#  Indoor Coverage with the PyLayers Multi Wall model


A Multi wall model accounts only for the attenuation along the direct path between Tx and Rx.


This approach do not provide information about delay spread or multi paths but it can nevertheless be useful in different contexts, as i.e optimization of an indoor radio network. The MultiWall approach provides a fast rough indication about the propagation losses.

A ray tracing approach is much more accurate, but  also is much more time consuming and depending on the purpose, it is relevant
 to proceed with a simpler and faster site-specific approach as the Multiwall model.


**`PyLayers`** provides a multiwall module which heavily relies on the core class **`Slab`**. Notice that, the same core **`Slab`** module is used for Ray tracing and MultiWall model approaches.

Let's see, how it works.
First let's import the `coverage` module. And the `time` module for performance evaluation.

```python
>>> from pylayers.antprop.coverage import *
>>> import time
```

Instantiate a coverage object.
By defaut the `TA-Office.ini` layout is loaded.

The coverage information is placed in the **coverage.ini** file in the project directory.

Below is an example such a configuration file.

```python
>>> !cat $BASENAME/ini/coverage.ini
[grid]
nx = 40
ny = 20
boundary = [20,0,30,20]
mode = full ; file, zone , full 
file = 'points.ini'

[layout]
filename = TA-Office.ini ; 0 40 0 15
;filename = W2PTIN.ini
;filename = Lstruc.str

[ap]
0 = {'name':'room1','wstd':'ieee80211b','p':(1,12,1.2),'PtdBm':0,'chan':[11],'on':True} 
1 = {'name':'room2','wstd':'ieee80211b','p':(10,2,1.2),'PtdBm':0,'chan':[11],'on':True} 
2 = {'name':'room3','wstd':'ieee80211b','p':(20,1,1.2),'PtdBm':0,'chan':[11],'on':True} 
3 = {'name':'room4','wstd':'ieee80211b','p':(36.5,1.5,1.2),'PtdBm':0,'chan':[11],'on':True} 
4 = {'name':'room5','wstd':'ieee80211b','p':(25,12,1.2),'PtdBm':0,'chan':[11],'on':True} 

[rx]
temperaturek = 300
noisefactordb = 0 

[show]
show = True
```

The ini file contains 5 sections.

+ [grid] section
    This section precises the size of the grid. By default the grid is placed over the whole region of the Layout.
    A selected region can also be defined whith the `boundary` list
+ [layout] section
    The name of the layout file (filename = )
+ [ap] section
    A dictionnary of access points precising the standard, the used channel, the emitted power and the position
    of the access point.

+ [show] section


```python
>>> # Create a Coverage object from coverag.ini file
... C = Coverage('coverage.ini')
```



```python
>>> C.dap[1]
name : room2
p : (10, 2, 1.2)
PtdBm : 0
channels  : [11]   2.462 : [2.451,2.473]
sensdBm : -94
nant : 1
On : True
```


The coverage object has a `__repr__` method which summarizes different parameters of the current coverage object

```python
>>> C
Layout file : TA-Office.ini

-----list of Access Points ------
name : room1
p : (1, 12, 1.2)
PtdBm : 0
channels  : [11]   2.462 : [2.451,2.473]
sensdBm : -94
nant : 1
On : True

name : room2
p : (10, 2, 1.2)
PtdBm : 0
channels  : [11]   2.462 : [2.451,2.473]
sensdBm : -94
nant : 1
On : True

name : room3
p : (20, 1, 1.2)
PtdBm : 0
channels  : [11]   2.462 : [2.451,2.473]
sensdBm : -94
nant : 1
On : True

name : room4
p : (36.5, 1.5, 1.2)
PtdBm : 0
channels  : [11]   2.462 : [2.451,2.473]
sensdBm : -94
nant : 1
On : True

name : room5
p : (25, 12, 1.2)
PtdBm : 0
channels  : [11]   2.462 : [2.451,2.473]
sensdBm : -94
nant : 1
On : True

-----Rx------
temperature (K) : 300
noisefactor (dB) : 0

--- Grid ----
mode : full
nx : 40
ny : 20
```

Then, the coverage calculation is launched by calling the `cover()` method


```python
>>> tic = time.time()
>>> C.cover()
>>> toc = time.time()
>>> print "Execution time : %2.3f " % (toc-tic)
Execution time : 1.687
```

Let display the current Layout with hidding nodes.

```python
>>> from matplotlib.pyplot import *
>>> %matplotlib inline
>>> fig=figure(figsize=(10,5))
>>> C.L.display['nodes']=False
>>> C.L.display['ednodes']=False
>>> f,a = C.show(fig=fig)
/home/uguen/anaconda/lib/python2.7/site-packages/matplotlib/collections.py:650: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors_original != str('face'):
/home/uguen/anaconda/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
```


The shadowing map coverage results can be displayed by invoquing various functions.

```python
>>> fig=figure(figsize=(10,5))
>>> f,a=C.show(fig=fig,typ='pr')
```

```python
>>> fig=figure(figsize=(10,5))
>>> f,a=C.show(fig=fig,typ='pr',f=4)
```

```python
>>> fig=figure(figsize=(10,5))
>>> f,a=C.show(fig=fig,typ='pr',f=10)
```

```python
>>> fig=figure(figsize=(10,5))
>>> f,a=C.show(fig=fig,typ='best',f=1)
```

```python
>>> fig=figure(figsize=(10,5))
>>> f,a=C.show(fig=fig,typ='best',f=10)
```

```python
>>> fig=figure(figsize=(10,5))
>>> C.show(fig=fig,f=5,typ='sinr')
(<matplotlib.figure.Figure at 0x7f1d80d48290>,
 <matplotlib.axes._subplots.AxesSubplot at 0x7f1d80b83850>)
```

As you have noticed the calculation has been done for all the center frequencies of the selected standard.
This is done in prevision of further channel optimizations.

Let's consider an other standard

```python
>>> C2 = Coverage('coverage2.ini')
>>> C2.cover()
```

```python
>>> fig=figure(figsize=(10,5))
>>> C2.show(ftyp='pr')
(<matplotlib.figure.Figure at 0x7f1d80fc2410>,
 <matplotlib.axes._subplots.AxesSubplot at 0x7f1d810945d0>)
```

```python
>>> C.snro.shape
(13, 800, 5)
```

```python
>>> fig=figure(figsize=(10,5))
>>> C.show(fig=fig,f=5,typ='capacity',dB=False)
(<matplotlib.figure.Figure at 0x7f1d8045c510>,
 <matplotlib.axes._subplots.AxesSubplot at 0x7f1d804541d0>)
```

All simulated quantities are stored in linear scale.

```python
>>> C2.Lwo[0,0,0]
0.078045027166146197
```

```python
>>> C2.freespace[0,0,0]
8.520773206944774e-07
```

