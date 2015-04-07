```python
>>> #!/usr/bin/python
... import time
>>> import matplotlib as ml
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> from types import *
>>> from numpy import array
>>> %pylab inline
>>> import pylayers.measures.vna.E5072A as E5072A
Populating the interactive namespace from numpy and matplotlib
```

```python
>>> vna = E5072A.SCPI("129.20.33.201",verbose=False)
```

```python
>>> # open remote measurement device (replace "hostname" by its actual name)
... data = vna.getIdent()
>>> print "Instrument ID : ",data
Instrument ID :
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
[<matplotlib.lines.Line2D at 0x7f283ef0f790>]
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
(50, 201)
```

```python
>>> R=np.mean(res,axis=0)
```

```python
>>> plt.plot(abs(R))
[<matplotlib.lines.Line2D at 0x7f283f45b590>]
```

```python
>>> r = ifft(R)
```

```python
>>> t = np.linspace(0,201/(2.2-1.8),201)
```

```python
>>> plt.plot(t*0.3,fftshift(abs(r)))
[<matplotlib.lines.Line2D at 0x7f283f05a190>]
```

```python
>>> plt.figure(figsize=(20,10))
>>> plt.imshow(abs(res),extent=(1.8,2.2,0,.1),origin='lower')
```

```python
>>> plt.plot(fftshift(abs(fres[0,:])))
[<matplotlib.lines.Line2D at 0x7f28444ff390>]
```

```python
>>> 3238-3216
22
```

```python
>>> len(S[22:])
3216
```

```python
>>> S21=np.frombuffer(S[0:201*16],dtype='>f8')
```

```python
>>> len(S21)
402
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
[<matplotlib.lines.Line2D at 0x7f2846de0950>]
```

```python
>>> plt.plot(freq,20*np.angle(cS21))
[<matplotlib.lines.Line2D at 0x7f2846ce23d0>]
```

```python
>>> import numpy as np
>>> f = np.frombuffer(tab,dtype='>i2')
```

```python
>>> 201*8
1608
```

```python
>>> fr=vna.getfreq()
```

```python
>>> S=vna.getnpoints()
```

```python
>>> vna.s.send(":SENS1:SWE:POIN?\n")
17
```

```python
>>> vna.s.recv(56)
'000000E+009,+1.80200000000E+009,+1.80400000000E+009,+1.8'
```

```python
>>> S=vna.getdata()
201
```

```python
>>> import pylayers.measures.switch.ni_usb_6501 as sw
>>> switch = sw.get_adapter()
>>> if not switch:
...     raise Exception("No device found")
>>> switch.set_io_mode(0b11111111, 0b11111111, 0b00000000)
'\x00\x08\x01\x00\x00\x00\x00\x02'
```

```python
>>> switch.write_port(0,0b00000101)
'\x00\x08\x01\x00\x00\x00\x00\x02'
```

```python
>>> eval('0b100')
4
```

```python

```
