# Using the vna

```python
>>> from E5072A import *
>>> from pylayers.antprop.channel import *
>>> %pylab inline
```

```python
>>> vna = SCPI("129.20.33.201",verbose=False)
```

```python
>>> vna.write(":SENS1:SWE:POIN 201")
```

```python
>>> vna.write(":SENS1:BAND 10")
```

```python
>>> vna = SCPI("129.20.33.201",verbose=False)
>>> ident = vna.getIdent()
>>> print "Talking to : ",ident
>>> vna.write("FORM:DATA REAL")
>>> vna.select(param='S21',chan=1)
>>> vna.write(":SENS1:SWE:POIN 1201")
>>> vna.write("DISP:WIND1:TRAC1:Y:SCAL:AUTO")
>>> vna.s.send(":SENS1:SWE:POIN?\n")
>>> Npoints = eval(vna.s.recv(56).replace('\n',''))
>>> print "Npoints : ",Npoints
>>> # set fmin fmax
... vna.write(":SENS1:FREQ:STAR 1.8e9")
>>> vna.write(":SENS1:FREQ:STOP 2.2e9")
>>> 
>>> #get frequency range
... com = ":SENS1:FREQ:DATA?\n"
>>> 
>>> vna.write("TRIG:SING")
>>> 
>>> time.sleep(1)
>>> com1 = ":CALC1:DATA:SDAT?\n"
>>> #u = np.arange(0,Npoints)*2
... #v = np.arange(0,Npoints)*2+1
... N = 1
>>> fGHz = np.linspace(1.8,2.2,1201)
>>> y  = np.zeros(len(fGHz))[None,:]
>>> H = Tchannel(x=fGHz,y=y)
>>> for k in range(N):
>>> #plt.imshow(abs(res))
...     B = vna.read(com1)
...     C = B[8:]
...     print len(C)
...     print Npoints*16
...     assert(len(C)==Npoints*16)
...     S = np.frombuffer(C,dtype='>f8')
...     H.frombuf(S)
...     #S21 = S[u]+1j*S[v]
...     #try:
...     #    res=np.vstack((res,S21.T))
...     #except:
...     #    res=S21.T
>>> #
... #
... #    tab = vna.readeval(vna)(com)
... #    f = np.frombuffer(tab,'>f8')
... #    freq = f[1:]
... #    plt.plot(freq)
... vna.close()
>>> #plt.imshow(abs(res))
... #plt.axis('tight')
... #plt.show()
```

```python
>>> H.save('calibration')
```

```python
>>> H.calibrate('calibration')
```

```python
>>> S
```

```python
>>> plt.plot(np.real(H.y[0]),'b',label="Real part")
>>> plt.plot(np.imag(H.y[0]),'r',label='Imaginary part"')
>>> plt.legend()
>>> plt.ylim(-1,2)
```

```python
>>> H.plot()
```

```python
>>> H
freq : 1.8 2.2 1201
shape  : (1, 1201)
tau (min, max) : [] []
dist :[] []

 calibrated : Yes
 windowed : No
```

```python
>>> fig=figure(figsize=(20,10))
>>> H.plot(fig=fig)
```

```python
>>> 800*0.1
80.0
```

```python
>>> h=H.ift(ffts=1)
```

```python
>>> h.plot(typ='v',xmin=-10,xmax=10)
```

```python
>>> H=Tchannel()
```

```python
>>> H.load('calibration.mat')
```

```python
>>> H
freq : 1.8 2.2 1201
shape  : (1, 1201)
tau (min, max) : [] []
dist :[] []

 calibrated : No
 windowed : No
```

```python
>>> H.plot()
(<matplotlib.figure.Figure at 0x7f9271740e90>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f9271737790>]], dtype=object))
```

```python
>>> H.calibrate('calibration.mat')
```

```python
>>> H.plot()
```

```python
>>> h=H.ift(ffts=1)
```

```python
>>> h.plot(typ='v')
```


