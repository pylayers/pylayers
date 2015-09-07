
Handling time and frequency domain signals : ``Bsignal`` Class
==============================================================

This section presents some features of the classes implemented in the
```pylayers.signal.bsignal.py`` <http://pylayers.github.io/pylayers/modules/pylayers.signal.bsignal.html>`__
module.

.. code:: python

    %matplotlib inline

The ``Bsignal`` class is a container for a signal with a base which can
be either in time domain or frequency domain.

.. code:: python

    from pylayers.signal.bsignal import *
    from matplotlib.pyplot import *

As a first example, let construct an impulse signal normalized in
energy. To do so there exist a specialized function :
```EnImpulse`` <http://pylayers.github.io/pylayers/modules/generated/pylayers.signal.bsignal.EnImpulse.demo.html#pylayers.signal.bsignal.EnImpulse.demo>`__

.. code:: python

    E=EnImpulse(fe=40)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-3-aa7406c5738b> in <module>()
    ----> 1 E=EnImpulse(fe=40)
    

    NameError: name 'EnImpulse' is not defined


.. code:: python

    print E.energy()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-4-8dbb63703b88> in <module>()
    ----> 1 print E.energy()
    

    NameError: name 'E' is not defined


.. code:: python

    E.plot(typ='v')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-5-9f7a16f65703> in <module>()
    ----> 1 E.plot(typ='v')
    

    NameError: name 'E' is not defined


.. code:: python

    E.energy()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-6-0fd52b52ad32> in <module>()
    ----> 1 E.energy()
    

    NameError: name 'E' is not defined


The Fourier transform of this signal has the hermitian Symmetry.

.. code:: python

    F = E.fft()
    F.plot(typ='m')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-7-d6500b4ff58a> in <module>()
    ----> 1 F = E.fft()
          2 F.plot(typ='m')


    NameError: name 'E' is not defined


.. code:: python

    F.y[0]


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-8-a3ec08f56226> in <module>()
    ----> 1 F.y[0]
    

    NameError: name 'F' is not defined


We then extract the non redundant part of the signal with the ``ft``
method

.. code:: python

    G=E.ft()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-9-d4e4ca97ebda> in <module>()
    ----> 1 G=E.ft()
    

    NameError: name 'E' is not defined


.. code:: python

    GH=G.symHz(100,scale='extract')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-10-e9dc0e1ba8d8> in <module>()
    ----> 1 GH=G.symHz(100,scale='extract')
    

    NameError: name 'G' is not defined


.. code:: python

    print GH.y[1]
    print GH.y[-1]


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-11-93f5cb675c6d> in <module>()
    ----> 1 print GH.y[1]
          2 print GH.y[-1]


    NameError: name 'GH' is not defined


.. code:: python

    ip=F.ifft()
    ip2=GH.ifft()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-12-a11b54673c68> in <module>()
    ----> 1 ip=F.ifft()
          2 ip2=GH.ifft()


    NameError: name 'F' is not defined


.. code:: python

    f,a=E.plot(typ='v',labels=['original'])
    f,a=ip.plot(typ='v',fig=f,ax=a[0][0],labels=['no zero padding'])
    f,a=ip2.plot(typ='v',fig=f,ax=a[0][0],labels=['zero padding'])
    title('extract mode')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-13-883c18d72357> in <module>()
    ----> 1 f,a=E.plot(typ='v',labels=['original'])
          2 f,a=ip.plot(typ='v',fig=f,ax=a[0][0],labels=['no zero padding'])
          3 f,a=ip2.plot(typ='v',fig=f,ax=a[0][0],labels=['zero padding'])
          4 title('extract mode')


    NameError: name 'E' is not defined


.. code:: python

    ip.energy()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-14-381d0b727c74> in <module>()
    ----> 1 ip.energy()
    

    NameError: name 'ip' is not defined


.. code:: python

    ip2.energy()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-15-78dac4d4a6c8> in <module>()
    ----> 1 ip2.energy()
    

    NameError: name 'ip2' is not defined


.. code:: python

    Y=E.esd()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-16-67b72beeca17> in <module>()
    ----> 1 Y=E.esd()
    

    NameError: name 'E' is not defined


FHsignal for in CIR mode
------------------------

We create a Fusignal which corresponds to the signal

.. math:: X_u(f) = \alpha e^{-2j\pi f \tau}

.. math:: f\in [f_{min},f_{max}]

.. code:: python

    f = np.arange(2,10,0.01)
    y = 2*np.ones(len(f))*np.exp(-2*1j*np.pi*f*3)
    N = len(f)
    Hu = FUsignal(f,y)
    print N


.. parsed-literal::

    800


.. code:: python

    Hu.plot(typ='m')




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2b78e57acdd0>,
     array([[<matplotlib.axes._subplots.AxesSubplot object at 0x2b78e580cc10>]], dtype=object))




.. image:: Bsignal_files/Bsignal_27_1.png


.. code:: python

    hu = Hu.ifft()

The inverse Fourier transform allows to recover perfectly the amplitude
:math:`\alpha` and the delay :math:`\tau` of the channel

.. code:: python

    hu.plot(typ='m')




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2b78e58e3650>,
     array([[<matplotlib.axes._subplots.AxesSubplot object at 0x2b78e5862850>]], dtype=object))




.. image:: Bsignal_files/Bsignal_30_1.png


.. code:: python

    real=np.imag(hu.y)
    u = np.where(hu.y==max(hu.y))[0]
    tau = hu.x[u]
    alpha = abs(hu.y[u])
    print alpha,tau


.. parsed-literal::

    [[  4.27928449e-14   4.45602724e-14   4.66296795e-14 ...,   3.83028329e-14
        3.93775414e-14   4.10699754e-14]
     [  4.27928449e-14   4.45602724e-14   4.66296795e-14 ...,   3.83028329e-14
        3.93775414e-14   4.10699754e-14]
     [  4.27928449e-14   4.45602724e-14   4.66296795e-14 ...,   3.83028329e-14
        3.93775414e-14   4.10699754e-14]
     ..., 
     [  4.27928449e-14   4.45602724e-14   4.66296795e-14 ...,   3.83028329e-14
        3.93775414e-14   4.10699754e-14]
     [  4.27928449e-14   4.45602724e-14   4.66296795e-14 ...,   3.83028329e-14
        3.93775414e-14   4.10699754e-14]
     [  4.27928449e-14   4.45602724e-14   4.66296795e-14 ...,   3.83028329e-14
        3.93775414e-14   4.10699754e-14]] [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.]


.. code:: python

    H = Hu.symHz(100,scale='cir')

.. code:: python

    H.plot(typ='m')




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2b78e59613d0>,
     array([[<matplotlib.axes._subplots.AxesSubplot object at 0x2b78e5ad04d0>]], dtype=object))




.. image:: Bsignal_files/Bsignal_33_1.png


.. code:: python

    h = H.ifft()

.. code:: python

    h.plot(typ='v')




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2b78e5a34d10>,
     array([[<matplotlib.axes._subplots.AxesSubplot object at 0x2b78e5a8ac90>]], dtype=object))




.. image:: Bsignal_files/Bsignal_35_1.png


.. code:: python

    real=np.imag(h.y)
    u = np.where(h.y==max(h.y))[0]
    tau = h.x[u]
    alpha = abs(h.y[u])
    print alpha,tau


.. parsed-literal::

    [[  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
        1.91361122e-02   1.16994371e-03]
     [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
        1.91361122e-02   1.16994371e-03]
     [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
        1.91361122e-02   1.16994371e-03]
     ..., 
     [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
        1.91361122e-02   1.16994371e-03]
     [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
        1.91361122e-02   1.16994371e-03]
     [  4.26036983e-14   5.77312005e-03   1.92248178e-02 ...,   2.34781108e-03
        1.91361122e-02   1.16994371e-03]] [-49.97728305 -49.97728305 -49.97728305 ..., -49.97728305 -49.97728305
     -49.97728305]


.. code:: python

    fft.ifft(H.y)




.. parsed-literal::

    array([[ -1.93565190e-15 -1.70240923e-19j,
              2.62295322e-04 -3.27871407e-19j,
              8.73458329e-04 -4.09839258e-20j, ...,
             -1.06670199e-04 +2.90350482e-19j,
             -8.69428086e-04 -1.58117458e-18j,
             -5.31550980e-05 -2.71727936e-20j]])



.. code:: python

    print H.y[203]
    print H.y[-203]
    len(H.y)


::


    ---------------------------------------------------------------------------

    IndexError                                Traceback (most recent call last)

    <ipython-input-28-4b5e7dcf9b6a> in <module>()
    ----> 1 print H.y[203]
          2 print H.y[-203]
          3 len(H.y)


    IndexError: index 203 is out of bounds for axis 0 with size 1


.. code:: python

    Y=h.fft()

.. code:: python

    Y.plot(typ='m')




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2b78e5c57ad0>,
     array([[<matplotlib.axes._subplots.AxesSubplot object at 0x2b78e5c57110>]], dtype=object))




.. image:: Bsignal_files/Bsignal_40_1.png

