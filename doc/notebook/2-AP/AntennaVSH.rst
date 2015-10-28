
Vector Spherical Harmonics Representation of Antennas
=====================================================

.. code:: python

    from pylayers.antprop.antenna import *
    from pylayers.antprop.antvsh import *
    %matplotlib inline


.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


Loading an Antenna from a Matlab file

.. code:: python

    A = Antenna('S2R2.mat',directory='ant/UWBAN/Matfile')

The shape of the :math:`F_{\phi}` functions indicates :

-  :math:`N_f= 104`
-  :math:`N_{\theta} = 91`
-  :math:`N_{\phi} = 180 `

.. code:: python

    np.shape(A.Fphi)


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-3-706142716bed> in <module>()
    ----> 1 np.shape(A.Fphi)
    

    AttributeError: 'Antenna' object has no attribute 'Fphi'


The frequency array is expressed in :math:`GHz` and delays are expressed
in :math:`ns`

.. code:: python

    fGHz = A.fa.reshape(104,1,1)


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-4-9c2d4af17123> in <module>()
    ----> 1 fGHz = A.fa.reshape(104,1,1)
    

    AttributeError: 'Antenna' object has no attribute 'fa'


Then an electrical delay of :math:`4.185ns` is applied on the
:math:`F_{\theta}`

.. code:: python

    I = A.Ftheta[:,:,:]
    plt.figure(figsize=(10,8))
    plt.imshow(np.unwrap(np.angle(I[:,45,:])))
    plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
    plt.ylabel('f index')
    plt.colorbar()
    plt.figure()
    plt.plot(fGHz[:,0,0],np.unwrap(np.angle(I[:,45,85])))
    plt.xlabel('f index')


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-5-087f4a8fccdb> in <module>()
    ----> 1 I = A.Ftheta[:,:,:]
          2 plt.figure(figsize=(10,8))
          3 plt.imshow(np.unwrap(np.angle(I[:,45,:])))
          4 plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
          5 plt.ylabel('f index')


    AttributeError: 'Antenna' object has no attribute 'Ftheta'


.. code:: python

    tau=4.185
    I = A.Ftheta[:,:,:]*exp(-2*1j*pi*fGHz*tau)


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-6-6b2a5143fb58> in <module>()
          1 tau=4.185
    ----> 2 I = A.Ftheta[:,:,:]*exp(-2*1j*pi*fGHz*tau)
    

    AttributeError: 'Antenna' object has no attribute 'Ftheta'


.. code:: python

    plt.imshow(np.unwrap(np.angle(I[:,45,:])))
    plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
    plt.ylabel('f index')
    plt.colorbar()
    
    plt.plot(fGHz[:,0,0],np.unwrap(np.angle(I[:,45,85])))


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-7-f05e7b6eaa51> in <module>()
    ----> 1 plt.imshow(np.unwrap(np.angle(I[:,45,:])))
          2 plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
          3 plt.ylabel('f index')
          4 plt.colorbar()
          5 


    NameError: name 'I' is not defined


Display of the radiation pattern for all frequencies
''''''''''''''''''''''''''''''''''''''''''''''''''''

.. code:: python

    plt.figure(figsize=(10,10))
    for nf in range(104):
        plt.polar(A.phi,abs(A.Ftheta[nf,45,:]))


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-8-dcc0dc965307> in <module>()
          1 plt.figure(figsize=(10,10))
          2 for nf in range(104):
    ----> 3     plt.polar(A.phi,abs(A.Ftheta[nf,45,:]))
    

    AttributeError: 'Antenna' object has no attribute 'Ftheta'



.. parsed-literal::

    <matplotlib.figure.Figure at 0x2af56d0b7510>


.. code:: python

    print 'Ntheta',A.Nt
    print 'Nphi',A.Np
    print 'Nf',A.Nf


.. parsed-literal::

    Ntheta

::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-9-9a627edc728e> in <module>()
    ----> 1 print 'Ntheta',A.Nt
          2 print 'Nphi',A.Np
          3 print 'Nf',A.Nf


    AttributeError: 'Antenna' object has no attribute 'Nt'


.. code:: python

    A.info()


.. parsed-literal::

     S2R2.mat
    type :  mat
    S2R2
    Th1
    04/13/12
    09:59
    
    
    2
    2
    Nb theta (lat) : 91
    Nb phi (lon) : 180
    No vsh coefficient calculated yet


Evaluation of Vector Spherical Harmonics Coefficients
=====================================================

At that stage we compute the Vector Spherical Harmonics coefficients

.. code:: python

    A=vsh(A)

.. code:: python

    A.info()


.. parsed-literal::

    S2R2.mat
    type :  mat
    S2R2
    Th1
    04/13/12
    09:59
    
    
    2
    2
    Nb theta (lat) : 91
    Nb phi (lon) : 180
    No vsh coefficient calculated yet


.. code:: python

    A.C.s1tos2(30)

.. code:: python

    A.C




.. parsed-literal::

    Br
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    
    Bi
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    
    Cr
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    
    Ci
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495



.. code:: python

    fig = plt.figure(figsize=(8,8))
    A.C.show('s2',k=300)



.. image:: AntennaVSH_files/AntennaVSH_23_0.png


.. code:: python

    A.C.s2tos3()

.. code:: python

    A.C




.. parsed-literal::

    Br
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    Ncoeff s3 : 145
    
    Bi
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    Ncoeff s3 : 145
    
    Cr
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    Ncoeff s3 : 145
    
    Ci
    -------------
    L1  : 90
    M1  : 89
    Ncoeff s1 8010
    NCoeff s2  : 495
    Ncoeff s3 : 145



.. code:: python

    fig = plt.figure(figsize=(8,8))
    A.C.show('s3')
    plt.tight_layout()



.. image:: AntennaVSH_files/AntennaVSH_26_0.png


