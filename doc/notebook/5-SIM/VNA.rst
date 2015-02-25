
.. code:: python

    #!/usr/bin/python
    import time
    import matplotlib as ml
    import matplotlib.pyplot as plt
    import numpy as np 
    from types import *
    from numpy import array
    %matplotlib inline
    import E5072A

::


    ---------------------------------------------------------------------------
    ImportError                               Traceback (most recent call last)

    <ipython-input-1-7b5128bbff9a> in <module>()
          7 from numpy import array
          8 get_ipython().magic(u'matplotlib inline')
    ----> 9 import E5072A
    

    ImportError: No module named E5072A


.. code:: python

    vna = E5072A.SCPI("129.20.33.201",verbose=False)
.. code:: python

    # open remote measurement device (replace "hostname" by its actual name)
    data = vna.getIdent()
    print "Instrument ID : ",data
.. code:: python

    vna.s.send('')
.. code:: python

    vna.select(param='S21',chan=1)
.. code:: python

    com =":SWE:POIN 1201"
    vna.write(com)
.. code:: python

    com = ":SENS1:FREQ:DATA?\n"
    tab = vna.read(com)
    f = np.frombuffer(tab,'>f8')
    freq = f[1:]
    plt.plot(freq)
.. code:: python

    try:
        del res
    except:
        pass
    com1 = "FORM:DATA REAL"
    com2 = "TRIG:SING"
    vna.write(com1)
    vna.write(com2)
    u = np.arange(0,201)*2
    v = np.arange(0,201)*2+1
    com = ":CALC1:DATA:SDAT?\n"
    N = 50
    for k in range(N):  
        B = vna.read(com)
        S =np.frombuffer(B[0:201*16],dtype='>f8')
        S21= S[u]+1j*S[v]
        try:
            res = np.vstack((res,S21.T))
        except:
            res = S21.T
.. code:: python

    from scipy.fftpack import fft,ifft,fftshift
.. code:: python

    fres=ifft(res,axis=1)
.. code:: python

    np.shape(res)
.. code:: python

    R=np.mean(res,axis=0)
.. code:: python

    plt.plot(abs(R))
.. code:: python

    r = ifft(R)
.. code:: python

    t = np.linspace(0,201/(2.2-1.8),201)
.. code:: python

    plt.plot(t*0.3,fftshift(abs(r)))
.. code:: python

    plt.figure(figsize=(20,10))
    plt.imshow(abs(res),extent=(1.8,2.2,0,.1),origin='lower')
.. code:: python

    plt.plot(fftshift(abs(fres[0,:])))
.. code:: python

    3238-3216
.. code:: python

    len(S[22:])
.. code:: python

    S21=np.frombuffer(S[0:201*16],dtype='>f8')
.. code:: python

    len(S21)
.. code:: python

    u = np.arange(0,201)*2
    v = np.arange(0,201)*2+1
.. code:: python

    cS21= S21[u]+1j*S21[v]
.. code:: python

    plt.plot(freq,20*np.log10(abs(cS21)))
.. code:: python

    plt.plot(freq,20*np.angle(cS21))
.. code:: python

    import numpy as np
    f = np.frombuffer(tab,dtype='>i2')
.. code:: python

    201*8
.. code:: python

    fr=vna.getfreq()
.. code:: python

    S=vna.getnpoints()
.. code:: python

    vna.s.send(":SENS1:SWE:POIN?\n")
.. code:: python

    vna.s.recv(56)
.. code:: python

    S=vna.getdata()
.. code:: python

    import pylayers.measures.switch.ni_usb_6501 as sw
    switch = sw.get_adapter()
    if not switch:
        raise Exception("No device found")
    switch.set_io_mode(0b11111111, 0b11111111, 0b00000000)

.. code:: python

    switch.write_port(0,0b00000101)
.. code:: python

    eval('0b100')
