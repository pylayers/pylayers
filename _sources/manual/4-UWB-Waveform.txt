.. code:: python

    from pylayers.signal.bsignal import *
    from pylayers.simul.simulem import *
    figsize(8,8)

Generation of an Impulse of normalized energy
=============================================

One possible manner to define an ernergy normalized short UWB impulse is
as follows.

:math:`p(t)= \frac{\sqrt{2\sqrt{2}}}{\tau\sqrt{\pi}} \cos(2\pi f_c t) e^{-(\frac{t}{\tau})^2}`

:math:`\tau = \frac{2}{B\pi}\sqrt{\frac{\gamma_{dB}\ln{10}}{20}}`

where :math:`B` is the desired bandwidth defined at
:math:`\gamma_{dB}` below the spectrum maximum and :math:`f_c` is
the central frequency of the pulse.

This waveform is a gaussian windowing of a sine wave of frequency
:math:`f_c`. The normalization term depends on the exponential scaling
factor :math:`\tau`.

.. code:: python

    fc     = 4 
    band   = 2
    thresh = 10
    fe     = 100 
    ip =EnImpulse([],fc,band,thresh,fe)

.. code:: python

    ip.info()

.. parsed-literal::

    TUsignal
    --------
    shx :  (343,)
    shy :  (343,)
    dx :   0.01
    xmin : -1.71
    xmax : 1.71
    ymin : -1.89545539648
    ymax : 2.16154131873


Verification of energy normalization in both domains
----------------------------------------------------


.. code:: python

    E1= sum(ip.y*ip.y)*ip.dx()
    print "Integration in time",E1

.. parsed-literal::

    Integration in time 1.0


.. code:: python

    P = ip.esd()
    E2 = sum(P.y)*P.dx()
    print "Integration in frequency domain ",E2

.. parsed-literal::

    Integration in frequency domain  1.0


Calcul of UWB channel impulse response
--------------------------------------



.. code:: python

    S= Simul()
    S.load('where2.ini')


.. code:: python

    st = S.wav.st
    sf = S.wav.sf
    S.wav.info()

.. parsed-literal::

    tw  :  30.0
    band  :  4.0
    fc  :  4.493
    thresh  :  3.0
    fe  :  50.0
    Np  :  1500.0
    te  :  0.02
    type  :  file


Here the time domain waveform is measured and the anticausal part of the
signal is artificially set to 0.

To handle properly the time domain wavefom it is required to center the
signal in the middle of the array.

``st`` stands for signal in time domain

Ploting the waveform
--------------------

In time domain
~~~~~~~~~~~~~~~

In[36]:

.. code:: python

    S.wav.st.plot()


.. image:: 4__UWB_Waveform_files/4__UWB_Waveform_fig_00.png

in frequency domain 
~~~~~~~~~~~~~~~~~~~

The frequency domain version of the signal is embedded in the same
object.

``sf`` stands for signal in frequency domain.


.. code:: python

    f,ax=S.wav.sf.plot()

.. image:: 4__UWB_Waveform_files/4__UWB_Waveform_fig_01.png

Construction of the propagation channel 
----------------------------------------

The link between Txid = 1 and Rxid =1 is simply loaded as


.. code:: python

    vc = S.VC(1,1)

.. parsed-literal::

    nray :  500
    nfreq :  181
    nb rays in .tauk file:  500
    nb rays 2:  500


The following representation shows the spatial spreading of the
propagation channel. On the left are scattered the intensity of rays wrt
to angles of departure (in azimut and elevation). On the right is the
intensity of rays wrt to angles of arrival. It misses the application
between the 2 planes as well as the delay dimension of the propagation
channel.


.. code:: python

    vc.doadod()

.. image:: 4__UWB_Waveform_files/4__UWB_Waveform_fig_02.png

Construction of the transmission channel
----------------------------------------


The transmission channel is obtain from the combianation of the
propagation channel and the vector antenna pattern at bot side of the
radio link


.. code:: python

    sc = vc.vec2scal()

The ScalChannel object contains all the information about the ray
transfer functions. The transmission channel is obtained by applying a
vector radiation pattern using an antenna file. In the presented case,
it comes from a real antenna which has been used during the FP7 WHERE1
measurement campaign M1.

.. code:: python

    S.tx.A.info()

.. parsed-literal::

    defant.vsh3
    type :  vsh3
    --------------------------
    fmin (GHz) : 2.0
    fmax (GHz) : 8.0
    Nf   : 121
    Br
    -------------
    Nf   :  121
    fmin (GHz) :  2.0
    fmax (GHz) :  8.0
    Ncoeff s3 :  18
    Bi
    -------------
    Nf   :  121
    fmin (GHz) :  2.0
    fmax (GHz) :  8.0
    Ncoeff s3 :  18
    Cr
    -------------
    Nf   :  121
    fmin (GHz) :  2.0
    fmax (GHz) :  8.0
    Ncoeff s3 :  18
    Ci
    -------------
    Nf   :  121
    fmin (GHz) :  2.0
    fmax (GHz) :  8.0
    Ncoeff s3 :  18



.. code:: python

    f,ax=sc.H.plot()


.. image:: 4__UWB_Waveform_files/4__UWB_Waveform_fig_03.png

The antenna can also been taken into account


.. code:: python

    alpha = 1./sqrt(30)  # scaling constant depends on how are stored the antenna data
    sca = vc.vec2scalA(S.tx.A,S.rx.A,alpha)
    sca.H.plot()

.. parsed-literal::

    (<matplotlib.figure.Figure at 0xb1c9b2c>,
     array([Axes(0.125,0.547727;0.775x0.352273),
           Axes(0.125,0.125;0.775x0.352273)], dtype=object))

.. image:: 4__UWB_Waveform_files/4__UWB_Waveform_fig_04.png

Calculate UWB Channel Impulse Response
--------------------------------------

.. code:: python

    cir = sc.applywavB(S.wav.sfg)

.. code:: python

    cir.plot()

.. image:: 4__UWB_Waveform_files/4__UWB_Waveform_fig_05.png

