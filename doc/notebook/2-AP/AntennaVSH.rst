
Vector Spherical Harmonics Representation of Antennas
=====================================================

.. code:: python

    >>> from pylayers.antprop.antenna import *
    >>> from pylayers.antprop.antvsh import *
    >>> %matplotlib inline

Loading an Antenna from a Matlab file

.. code:: python

    A = Antenna('S2R2.mat',directory='ant/UWBAN/Matfile')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-1-806a96154428> in <module>()
    ----> 1 A = Antenna('S2R2.mat',directory='ant/UWBAN/Matfile')
    

    NameError: name 'Antenna' is not defined


The shape of the :math:`F_{\phi}` functions indicates :

-  :math:`N_f= 104`
-  :math:`N_{\theta} = 91`
-  :math:`N_{\phi} = 180 `

.. code:: python

    np.shape(A.Fp)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-2-1104363ae4f5> in <module>()
    ----> 1 np.shape(A.Fp)
    

    NameError: name 'np' is not defined


The frequency array is expressed in :math:`GHz` and delays are expressed
in :math:`ns`

.. code:: python

    fGHz = A.fGHz


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-3-7837d8426b1a> in <module>()
    ----> 1 fGHz = A.fGHz
    

    NameError: name 'A' is not defined


.. code:: python

    fGHz.shape


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-4-9aefd7d028df> in <module>()
    ----> 1 fGHz.shape
    

    NameError: name 'fGHz' is not defined


Then an electrical delay of :math:`4.185ns` is applied on the
:math:`F_{\theta}`

.. code:: python

    I = A.Ft[:,:,:]


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-5-c10f1b8c978a> in <module>()
    ----> 1 I = A.Ft[:,:,:]
    

    NameError: name 'A' is not defined


.. code:: python

    I.shape


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-6-d1198b51cdae> in <module>()
    ----> 1 I.shape
    

    NameError: name 'I' is not defined


.. code:: python

    
    plt.figure(figsize=(10,8))
    plt.imshow(np.unwrap(np.angle(I[:,45,:])))
    plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
    plt.ylabel('f index')
    plt.colorbar()
    plt.figure()
    plt.plot(fGHz,np.unwrap(np.angle(I[45,85,:])))
    plt.xlabel('f index')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-7-790e88a91dc8> in <module>()
          1 
    ----> 2 plt.figure(figsize=(10,8))
          3 plt.imshow(np.unwrap(np.angle(I[:,45,:])))
          4 plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
          5 plt.ylabel('f index')


    NameError: name 'plt' is not defined


.. code:: python

    tau=4.185
    I = A.Ft[:,:,:]*np.exp(-2*1j*np.pi*fGHz[None,None,:]*tau)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-8-24216f96e669> in <module>()
          1 tau=4.185
    ----> 2 I = A.Ft[:,:,:]*np.exp(-2*1j*np.pi*fGHz[None,None,:]*tau)
    

    NameError: name 'A' is not defined


.. code:: python

    plt.imshow(np.unwrap(np.angle(I[:,45,:])))
    plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
    plt.ylabel('f index')
    plt.colorbar()
    plt.figure()
    plt.plot(fGHz,np.unwrap(np.angle(I[45,85,:])))


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-9-510d417ac34a> in <module>()
    ----> 1 plt.imshow(np.unwrap(np.angle(I[:,45,:])))
          2 plt.title(r'Unwrapped phase of $F_{\theta}$ w.r.t frequency and phi for $\theta=\frac{pi}{2}$')
          3 plt.ylabel('f index')
          4 plt.colorbar()
          5 plt.figure()


    NameError: name 'plt' is not defined


Display of the radiation pattern for all frequencies
''''''''''''''''''''''''''''''''''''''''''''''''''''

.. code:: python

    plt.figure(figsize=(10,10))
    for nf in range(104):
        plt.polar(A.phi,abs(A.Ft[45,:,nf]))


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-10-f991f6d8f1a6> in <module>()
    ----> 1 plt.figure(figsize=(10,10))
          2 for nf in range(104):
          3     plt.polar(A.phi,abs(A.Ft[45,:,nf]))


    NameError: name 'plt' is not defined


.. code:: python

    A.info()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-11-a01e925d9d7e> in <module>()
    ----> 1 A.info()
    

    NameError: name 'A' is not defined


Evaluation of Vector Spherical Harmonics Coefficients
=====================================================

At that stage we compute the Vector Spherical Harmonics coefficients

.. code:: python

    A=vsh(A)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-12-aab26d118494> in <module>()
    ----> 1 A=vsh(A)
    

    NameError: name 'vsh' is not defined


.. code:: python

    A.info()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-13-a01e925d9d7e> in <module>()
    ----> 1 A.info()
    

    NameError: name 'A' is not defined


.. code:: python

    A.C.s1tos2(30)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-14-98aa5dfbfef3> in <module>()
    ----> 1 A.C.s1tos2(30)
    

    NameError: name 'A' is not defined


.. code:: python

    A.C


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-15-ea02b37ef526> in <module>()
    ----> 1 A.C
    

    NameError: name 'A' is not defined


.. code:: python

    fig = plt.figure(figsize=(8,8))
    A.C.show('s2',k=300)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-16-a2dac715dfe4> in <module>()
    ----> 1 fig = plt.figure(figsize=(8,8))
          2 A.C.show('s2',k=300)


    NameError: name 'plt' is not defined


.. code:: python

    A.C.s2tos3()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-17-34ddad199ddd> in <module>()
    ----> 1 A.C.s2tos3()
    

    NameError: name 'A' is not defined


.. code:: python

    A.C


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-18-ea02b37ef526> in <module>()
    ----> 1 A.C
    

    NameError: name 'A' is not defined


.. code:: python

    fig = plt.figure(figsize=(8,8))
    A.C.show('s3')
    plt.tight_layout()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-19-627adf1c1577> in <module>()
    ----> 1 fig = plt.figure(figsize=(8,8))
          2 A.C.show('s3')
          3 plt.tight_layout()


    NameError: name 'plt' is not defined


