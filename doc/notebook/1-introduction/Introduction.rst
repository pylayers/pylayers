
What PyLayers is about ?
------------------------
<<<<<<< HEAD


-  Indoor Radio Propagation Modelling
-  IR UWB (IEEE 802.15.4a and IEEE 802.15.6)
-  Wireless Body Area Network
-  Indoor Localisation
-  Human mobility simulation

-  Cooperative Positioning

=======
>>>>>>> 964dfcf95e3ef8fe011e741d5de1b64b65094db6

-  Indoor Radio Propagation Modelling
-  IR UWB (IEEE 802.15.4a and IEEE 802.15.6)
-  Wireless Body Area Network
-  Indoor Localisation
-  Human mobility simulation

-  Cooperative Positioning

How to install PyLayers ?
-------------------------

<<<<<<< HEAD

=======
>>>>>>> 964dfcf95e3ef8fe011e741d5de1b64b65094db6
If you are working on a Windows platform and you are not familiar with
the Python ecosystem, a good idea is to install first
`Anaconda <https://store.continuum.io/cshop/anaconda/>`_. It will
provide you most of the required dependencies.

Dependencies
~~~~~~~~~~~~

::

    Before to install pylayers you need to install the following dependencies 


    numpy>=1.6.1
    scipy>=0.10.1
    networkx>=1.7
    matplotlib>=1.1.0
    shapely>=1.2.14
    descartes>=1.0
    SimPy>=2.2
    PIL>=1.1.5
    bitstring>=3.0.2
    pyintervall>=1.0
    pandas >=0.12.0

    The installation of matplotlib basemap requires installation of libraries GEOS
    and PROJ4. The installation procedure is well described in

http://matplotlib.org/basemap/users/installing.html

::

    For supporting the osm format PyLayers is relying on imposm. The installation 
    of imposm is not straighforward it requires to install first Tokyo-Cabinet
    and Google-Protobuf library. PyLayers can work without imposm you need to 
    comment the only import to the module in pylayers/gis/osmparser.py. 

http://imposm.org/docs/imposm/latest/install.html#id1

::

    Most of these dependencies are handled in the setup.py script.

    python setup.py install

Setting of environment variables $BASENAME and $PYLAYERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    It is required to define the two enviroment variables $PYLAYERS and $BASENAME

    in your .bashrc add those 2 lines 

    #export PYLAYERS = <path to the directory where pylayers is installed>

    for example if you install PyLyers in your home directory 

    export PYLAYERS=~/pylayers
    export BASENAME=~/plproject

Installing PyInterval
~~~~~~~~~~~~~~~~~~~~~

::

    get crlibm library first 
    http://lipforge.ens-lyon.fr/frs/?group_id=8&release_id=123

    if problem for compilation add 
    CPPFLAGS= -fPIC in scs_lib/Makefile

    sudo easy_install pyinterval

Testing
~~~~~~~

make test-code
