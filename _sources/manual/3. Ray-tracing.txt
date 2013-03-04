.. warning::
    This demo works only on Linux platform
    platform. It relies on exec files (launching, tracing,tratotud,evalfield) in the `bin` 
    directory of the project. It is foreseen to replace soon those dependencies by
    pure python equivalents. 

.. code:: python

    from pylayers.simul.simulem import *

Ray-tracing
-----------


.. code:: python

    S = Simul()
    # loading a layout 
    filestr = 'TA-Office'
    # load the material and slab dictionnaries
    S.layout(filestr+'.str','matDB.ini','slabDB.ini')

Setting the transmitter
~~~~~~~~~~~~~~~~~~~~~~~


.. code:: python

    S.tx = RadioNode(typ='tx')
    S.tx.point([1.2,1,1.4])

Setting receiver
~~~~~~~~~~~~~~~~



.. code:: python

    S.rx = RadioNode(typ='rx')
    S.rx.point([8,-1.2,1.5],mode='append')
    S.rx.point([8,-1.21,1.5],mode='append')
    S.rx.point([8,-1.22,1.5],mode='append')
    S.rx.point([8,-1.23,1.5],mode='append')
    S.rx.point([8,-1.24,1.5],mode='append')
    
    S.save()

Adjust simulation parameter (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



.. code:: python

    
    # print launching parameters
    S.palch.info()
    
    # ang Tx : angular step from Tx
    S.palch.angTx  = 1
    
    # ISB ang Incident Shadow Boundary angle (degree) 
    S.palch.ISBang = 90  
    
    # ray elimination Threshold 
    S.palch.ethreshold = 0.001
    
    # maximum depth
    S.palch.maxdeep  = 10
    
    # typealgo = 0 (include diffraction) 1 (no diffraction)
    S.palch.typalgo = 1
    title = str(S.palch.angTx) + '-' +\
            str(S.palch.ISBang) + '-' +\
            str(S.palch.ethreshold) + '-' + \
            str(S.palch.maxdeep) + '-' + \
            str(S.palch.typalgo)
    
    S.palch.save()
    S.pafreq.fghzmin=2
    S.pafreq.fghzmax=11
    S.pafreq.nf=181
    S.pafreq.save()

.. parsed-literal::

    ----------------------------------------------
                Launching Parameter               
    ----------------------------------------------
    angTx      : Tx angular step ( degrees)     :  1
    ISBang     : ISB angular sector ( degrees ) :  90
    ethreshold : Exploration Threshold (linear) :  0.001
    maxdeep    : Tree deep max (integer value)  :  10
    typalgo    : Type of algo (default 0)       :  1


Run the ray tracing
~~~~~~~~~~~~~~~~~~~



.. code:: python

    S.run(1,1)

.. parsed-literal::

    run debug  1 1
    nray : 

.. parsed-literal::

     500
    nfreq :  181
    nb rayons dans .tauk : 

.. parsed-literal::

     500
    nb rayons 2:  500


.. parsed-literal::

    /usr/local/lib/python2.7/site-packages/scipy-0.11.0-py2.7-linux-i686.egg/scipy/io/matlab/mio.py:266: FutureWarning: Using oned_as default value ('column') This will change to 'row' in future versions
      oned_as=oned_as)



.. parsed-literal::

    True

Plot CIR
~~~~~~~~



.. code:: python

    S.pltcir(1,1)

.. parsed-literal::

    Warning : no furniture file loaded


.. image:: 3__Ray_tracing_files/3__Ray_tracing_fig_00.png
