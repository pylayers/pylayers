
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
.. code:: python

    >>> from pylayers.signal.bsignal import *
    >>> ip    = EnImpulse(fc=4,band=3,thresh=10,fe=100)
    >>> Eip1  = ip.energy()
    >>> ESDu  = ip.esd(mode='unilateral')
    >>> ESDb  = ip.esd(mode='bilateral')
    >>> df    = ESDu.dx()
    >>> Eipu  = sum(ESDu.y)*df
    >>> Eipb  = sum(ESDb.y)*df
    >>> erru  = Eip1-Eipu
    >>> errb  = Eip1-Eipb
.. code:: python

    print Eip1

.. parsed-literal::

    100.000007743


.. code:: python

    E.plot(typ='v')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2ac10ceedb50>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10cef8090>]], dtype=object))




.. image:: Bsignal_files/Bsignal_9_1.png


.. code:: python

    E.energy()



.. parsed-literal::

    40.000003097054716



The Fourier transform of this signal has the hermitian Symmetry.

.. code:: python

    F = E.fft()
    F.plot(typ='m')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2ac10cf5b310>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10cf51050>]], dtype=object))




.. image:: Bsignal_files/Bsignal_12_1.png


.. code:: python

    F.y[0]



.. parsed-literal::

    (0.00029728997658457938+0j)



We then extract the non redundant part of the signal with the ``ft``
method

.. code:: python

    G=E.ft()
.. code:: python

    GH=G.symHz(100,scale='extract')
.. code:: python

    print GH.y[1]
    print GH.y[-1]

.. parsed-literal::

    (-0.0014441784194-4.88037298122e-05j)
    (-0.0014441784194+4.88037298122e-05j)


.. code:: python

    ip=F.ifft()
    ip2=GH.ifft()
.. code:: python

    f,a=E.plot(typ='v',labels=['original'])
    f,a=ip.plot(typ='v',fig=f,ax=a[0][0],labels=['no zero padding'])
    f,a=ip2.plot(typ='v',fig=f,ax=a[0][0],labels=['zero padding'])
    title('extract mode')



.. parsed-literal::

    <matplotlib.text.Text at 0x2ac10d0b6210>




.. image:: Bsignal_files/Bsignal_19_1.png


.. code:: python

    ip.energy()



.. parsed-literal::

    40.000003097054773



.. code:: python

    ip2.energy()



.. parsed-literal::

    401.35111446263471



.. code:: python

    Y=E.esd()


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

    (<matplotlib.figure.Figure at 0x2ac10cf71910>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10d02d350>]], dtype=object))




.. image:: Bsignal_files/Bsignal_27_1.png


.. code:: python

    hu = Hu.ifft()
The inverse Fourier transform allows to recover perfectly the amplitude
:math:`\alpha` and the delay :math:`\tau` of the channel

.. code:: python

    hu.plot(typ='m')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2ac10d38b610>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10d29c1d0>]], dtype=object))




.. image:: Bsignal_files/Bsignal_30_1.png


.. code:: python

    real=np.imag(hu.y)
    u = np.where(hu.y==max(hu.y))[0]
    tau = hu.x[u]
    alpha = abs(hu.y[u])
    print alpha,tau

.. parsed-literal::

    [ 2.] [ 3.00375469]


.. code:: python

    H = Hu.symHz(100,scale='cir')
.. code:: python

    H.plot(typ='m')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2ac10d423d90>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10d3c01d0>]], dtype=object))




.. image:: Bsignal_files/Bsignal_33_1.png


.. code:: python

    h = H.ifft()
.. code:: python

    h.plot(typ='v')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2ac10d4ef9d0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10d4e4b10>]], dtype=object))




.. image:: Bsignal_files/Bsignal_35_1.png


.. code:: python

    real=np.imag(h.y)
    u = np.where(h.y==max(h.y))[0]
    tau = h.x[u]
    alpha = abs(h.y[u])
    print alpha,tau

.. parsed-literal::

    [ 1.97995425] [-46.97864607]


.. code:: python

    fft.ifft(H.y)



.. parsed-literal::

    array([ -1.93565190e-15 -1.70240923e-19j,
             2.62295322e-04 -3.27871407e-19j,
             8.73458329e-04 -4.09839258e-20j, ...,
            -1.06670199e-04 +2.90350482e-19j,
            -8.69428086e-04 -1.58117458e-18j,  -5.31550980e-05 -2.71727936e-20j])



.. code:: python

    print H.y[203]
    print H.y[-203]
    len(H.y)

.. parsed-literal::

    (0.116169256529-0.0459946624208j)
    (0.116169256529+0.0459946624208j)




.. parsed-literal::

    2201



.. code:: python

    Y=h.fft()
.. code:: python

    Y.plot(typ='m')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2ac10d5ac1d0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x2ac10d5a24d0>]], dtype=object))




.. image:: Bsignal_files/Bsignal_40_1.png

