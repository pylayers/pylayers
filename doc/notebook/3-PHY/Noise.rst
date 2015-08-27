
Modelisation of the Thermal Noise
=================================

.. code:: python

    from pylayers.signal.bsignal import *
    %matplotlib inline

The bsignal module has a dedicated class for handling noise signal. To
create a white noise just type :

.. code:: python

    w = Noise()

The representation of the noise object provides information about
default values. In digital representation of noise the sampling
frequency is important. The noise signal is generated from a time
:math:`t_i` to a time :math:`t_f = t_i+T`. The default power spectral
density is :math:`-174dBm/Hz` and can be modified with the argument
``PSDdBmpHz``.

.. code:: python

    w




.. parsed-literal::

    Sampling frequency : 50 GHz
    ti  : 0ns 
    tf  : 100ns 
    ts  : 0.02ns 
    N   : 5000
    -------------
    DSP : -174 dBm/Hz
        : 3.98107170553e-21 Joules
    -------------
    Noise Figure : 0 dB
    Vrms : 9.97631157484e-05 Volts
    Variance : 9.94012046425e-09 V^2
    Power (dBm) /50 Ohms : -157.010299957 dBm
    Power realized /50 Ohms : -157.015783567 dBm



.. code:: python

    f,a=w.plot(typ='v')



.. image:: Noise_files/Noise_6_0.png


.. code:: python

    w.psd()


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-5-638d26c14f9e> in <module>()
    ----> 1 w.psd()
    

    AttributeError: 'Noise' object has no attribute 'psd'


.. code:: python

    w2 = w.fgating(fcGHz=4,BGHz=3)


::


    ---------------------------------------------------------------------------

    IndexError                                Traceback (most recent call last)

    <ipython-input-6-0901d96562a6> in <module>()
    ----> 1 w2 = w.fgating(fcGHz=4,BGHz=3)
    

    /home/uguen/Documents/rch/devel/pylayers/pylayers/signal/bsignal.pyc in fgating(self, fcGHz, BGHz, window)
       3616         else:
       3617             parity = 1
    -> 3618         U = N.unrex()
       3619         f = U.x
       3620         f1 = fcGHz - BGHz / 2.


    /home/uguen/Documents/rch/devel/pylayers/pylayers/signal/bsignal.pyc in unrex(self)
       3411         if np.mod(N, 2) == 0:
       3412             xu = self.x[1:(N + 2) / 2]
    -> 3413             yu = self.y[:,1:(N + 2) / 2]
       3414         # odd case
       3415         else:


    IndexError: too many indices for array


.. code:: python

    W2=w2.psd()
    W2.plotdB(mask=True)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-7-00fd8114e3d7> in <module>()
    ----> 1 W2=w2.psd()
          2 W2.plotdB(mask=True)


    NameError: name 'w2' is not defined


.. code:: python

    w.plot(typ='v')




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x2b9921022f10>,
     array([[<matplotlib.axes._subplots.AxesSubplot object at 0x2b9920ede390>]], dtype=object))




.. image:: Noise_files/Noise_10_1.png


.. code:: python

    ip=EnImpulse(fc=4.4928,band=0.4992,fe=100)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-9-307076f57f86> in <module>()
    ----> 1 ip=EnImpulse(fc=4.4928,band=0.4992,fe=100)
    

    NameError: name 'EnImpulse' is not defined


.. code:: python

    fig = plt.figure(figsize=(10,10))
    for k,snr in enumerate(range(30,-30,-10)):
        a = fig.add_subplot(3,2,k+1)
        ipn,n=ip.awgn(snr=snr,typ='snr')
        ipn.plot(typ='v',fig=fig,ax=a)
        a.set_title('SNR :'+str(snr)+' dB')
    plt.tight_layout()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-10-897bc488bfef> in <module>()
          2 for k,snr in enumerate(range(30,-30,-10)):
          3     a = fig.add_subplot(3,2,k+1)
    ----> 4     ipn,n=ip.awgn(snr=snr,typ='snr')
          5     ipn.plot(typ='v',fig=fig,ax=a)
          6     a.set_title('SNR :'+str(snr)+' dB')


    NameError: name 'ip' is not defined



.. image:: Noise_files/Noise_12_1.png


