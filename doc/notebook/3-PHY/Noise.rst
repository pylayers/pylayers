
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
    -------------
    DSP : -174 dBm/Hz
    NF : 0 dB
    Vrms : 9.97631157484e-05 Volts
    Variance : 1.00178020455e-08 V^2
    Power /50 Ohms : -157.010299957 dBm
    Power realized /50 Ohms : -156.981975587 dBm



.. code:: python

    f,a=w.plot(typ='v')


.. image:: Noise_files/Noise_6_0.png


.. code:: python

    w.psd()



.. parsed-literal::

    FUsignal :  (2500,)  (2500,) 
    Frequency (GHz) : 2500



.. code:: python

    w2 = w.fgating(fcGHz=4,BGHz=3)
.. code:: python

    W2=w2.psd()
    W2.plotdB(mask=True)


.. image:: Noise_files/Noise_9_0.png


.. code:: python

    w.plot(typ='v')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f3fc45fd590>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f3fc45f5d10>]], dtype=object))




.. image:: Noise_files/Noise_10_1.png


.. code:: python

    ip=EnImpulse(fe=100)
.. code:: python

    fig = plt.figure(figsize=(10,10))
    for k,snr in enumerate(range(30,-30,-10)):
        a = fig.add_subplot(3,2,k+1)
        ipn=ip.awgn(snr=snr,typ='snr')
        ipn.plot(typ='v',fig=fig,ax=a)
        a.set_title('SNR :'+str(snr)+' dB')
    plt.tight_layout()


.. image:: Noise_files/Noise_12_0.png

