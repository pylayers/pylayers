```python
>>> #!/usr/bin/python
... import time
>>> import matplotlib as ml
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> from types import *
>>> from numpy import array
>>> %matplotlib inline
>>> import pylayers.measures.vna.E5072A
```

```python
>>> vna = E5072A.SCPI("129.20.33.201",verbose=False)
```

```python
>>> # open remote measurement device (replace "hostname" by its actual name)
... data = vna.getIdent()
>>> print "Instrument ID : ",data
```

```python
>>> vna.s.send('')
```

```python
>>> vna.select(param='S21',chan=1)
```

```python
>>> com =":SWE:POIN 1201"
>>> vna.write(com)
```

```python
>>> com = ":SENS1:FREQ:DATA?\n"
>>> tab = vna.read(com)
>>> f = np.frombuffer(tab,'>f8')
>>> freq = f[1:]
>>> plt.plot(freq)
```

```python
>>> try:
...     del res
>>> except:
...     pass
>>> com1 = "FORM:DATA REAL"
>>> com2 = "TRIG:SING"
>>> vna.write(com1)
>>> vna.write(com2)
>>> u = np.arange(0,201)*2
>>> v = np.arange(0,201)*2+1
>>> com = ":CALC1:DATA:SDAT?\n"
>>> N = 50
>>> for k in range(N):
...     B = vna.read(com)
...     S =np.frombuffer(B[0:201*16],dtype='>f8')
...     S21= S[u]+1j*S[v]
...     try:
...         res = np.vstack((res,S21.T))
...     except:
...         res = S21.T
```

```python
>>> from scipy.fftpack import fft,ifft,fftshift
```

```python
>>> fres=ifft(res,axis=1)
```

```python
>>> np.shape(res)
```

```python
>>> R=np.mean(res,axis=0)
```

```python
>>> plt.plot(abs(R))
```

```python
>>> r = ifft(R)
```

```python
>>> t = np.linspace(0,201/(2.2-1.8),201)
```

```python
>>> plt.plot(t*0.3,fftshift(abs(r)))
```

```python
>>> plt.figure(figsize=(20,10))
>>> plt.imshow(abs(res),extent=(1.8,2.2,0,.1),origin='lower')
```

```python
>>> plt.plot(fftshift(abs(fres[0,:])))
```

```python
>>> 3238-3216
```

```python
>>> len(S[22:])
```

```python
>>> S21=np.frombuffer(S[0:201*16],dtype='>f8')
```

```python
>>> len(S21)
```

```python
>>> u = np.arange(0,201)*2
>>> v = np.arange(0,201)*2+1
```

```python
>>> cS21= S21[u]+1j*S21[v]
```

```python
>>> plt.plot(freq,20*np.log10(abs(cS21)))
```

```python
>>> plt.plot(freq,20*np.angle(cS21))
```

```python
>>> import numpy as np
>>> f = np.frombuffer(tab,dtype='>i2')
```

```python
>>> 201*8
```

```python
>>> fr=vna.getfreq()
```

```python
>>> S=vna.getnpoints()
```

```python
>>> vna.s.send(":SENS1:SWE:POIN?\n")
```

```python
>>> vna.s.recv(56)
```

```python
>>> S=vna.getdata()
```

```python
>>> import pylayers.measures.switch.ni_usb_6501 as sw
>>> switch = sw.get_adapter()
>>> if not switch:
...     raise Exception("No device found")
>>> switch.set_io_mode(0b11111111, 0b11111111, 0b00000000)
```

```python
>>> switch.write_port(0,0b00000101)
```

```python
>>> eval('0b100')
```

```python

```
