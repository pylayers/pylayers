
.. code:: python

    import pylayers.gis.ezone as ez
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    %matplotlib inline
.. code:: python

    prefix=ez.enctile(-1.5,47.5)
    print prefix

.. parsed-literal::

    N47W002


.. code:: python

    E=ez.Ezone(prefix)
.. code:: python

    E



.. parsed-literal::

    N47W002
    --------
    (-2, -1, 47, 48)
    latlon : [ -0.000 75118.790 cartesian :0.000 111194.505 ]



.. code:: python

    E.getdem()

.. parsed-literal::

    Load srtm file
    no aster file for this point


.. code:: python

    E.show(clim=[-20,180])



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7ffa46a9de50>,
     <matplotlib.axes.AxesSubplot at 0x7ffa46ab62d0>)




.. image:: building_files/building_5_1.png


.. code:: python

    E



.. parsed-literal::

    N47W002
    --------
    (-2, -1, 47, 48)
    latlon : [ -0.000 75118.790 cartesian :0.000 111194.505 ]



