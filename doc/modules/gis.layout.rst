.. _bsignal:

layout
======

This class handles the description of the indoor propagation environment. 
It manages various input and output data format. 


The most basic description of a Layout is given in a `.ini` file 
which contains the points and edges list of the Layout.

This `.ini` file can be edited via the interactive command 

.. code-block:: python 
    L.editor()

Then different graphs are processed in order to enrich the data structure for 
further computations as the calculation of signatures.

.. code-block:: python 
    L.build()

As the building can be computed once there is a command for writing anf
reading the layout additional information from gpickle format files. 

gpickle files are stored in $BASENAME/struc/ directory

.. code-block:: python 
    try:
        L.dumpr()
    except:
        L.build()
        L.dumpw()

The `.ini` file can be converted in a `.str2` format which is used for the 
C implementation of yhe ray tracing. This format is going to be deprecated. 

The `.str2` file is processed to generate a `.str` file which 
contains additional visibility relationships between vertices and edges. 

.. autoclass:: pylayers.gis.layout
    :members:

.. automodule:: pylayers.gis.layout
    :members:
