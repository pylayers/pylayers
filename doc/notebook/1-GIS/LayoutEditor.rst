
The Layout class
----------------

Introduction
~~~~~~~~~~~~

This section explains the main features of the **Layout** class.

A Layout is a representation of a floorplan, it is handled by the module
**pylayers.gis.layout**.

This module recognizes different file formats including geo-referenced
files in open street map format ``.osm``.

Using osm allows to take advantage of a mature floor plan editor *JOSM*
with the plugin *PicLayer*. This is well described in
http://wiki.openstreetmap.org/wiki/IndoorOSM

The ``pylayers.gis.osmparser`` module parses ``osm`` files.

See the following methods of the layout object

-  ``loadosm()``
-  ``saveosm()``

Structure of a Layout
~~~~~~~~~~~~~~~~~~~~~

At first a Layout is described by a set of points (negative index) and a
set of segments (positive index).

Points and segments are nodes of the :math:`\mathcal{G}_s` graph.

It is required to respect a strict non overlapping rule. **No segments
can recover partially or totally an other segment**.

This rule allows to capture topological relations of the network which
are exploited for further analysis.

Subsegments
~~~~~~~~~~~

To describe doors and windows, the concept of ``subsegment`` is
introduced.

A ``segment`` has attributes :

-  ``name`` : slab name
-  ``z`` : tuple of minimum and maximum heights with respect to ground
   (meters)
-  ``transition`` : a boolean indicating if a human can cross this
   segment. For example, segments associated with a door are transition
   segments but we will see later that it may be judicious to split
   space with transparent segments which have the name 'AIR'. Those
   segments are also ``transition=True``

A ``subsegment`` belongs to a ``segment``, it has mainly 2 attached
parameters :

-  ``ss_name`` : subsegment slab name
-  ``ss_z`` : [(zmin1,zmax1),(zmin2,zmax2),...,(zminK,zmaxK))] list of
   minimum and maximum height of associated subsegments (meters)

When appearing in a 3D ray a subsegment should have a unique index
different from the segment index.

The different layout format
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The layout format has regularly evolved over time and is going to evolve
again. Currently, the different recognized file extensions are the
following :

-  ``.str2``: a ASCII file (Node list + edge list)
-  ``.str`` : a binary file which includes visibility relations between
   point and segments
-  ``.ini`` : an ini file which gather node list and edge list as well
   as the state of the current ``display`` dictionnary
-  ``.osm`` : an xml file which can be edited with
   `JOSM <http://josm.openstreetmap.de/>`__

.. code:: python

    from pylayers.gis.loyout import *
    from pylayers.util.project import *


::


    ---------------------------------------------------------------------------

    ImportError                               Traceback (most recent call last)

    <ipython-input-1-341e1dead720> in <module>()
    ----> 1 from pylayers.gis.loyout import *
          2 from pylayers.util.project import *


    ImportError: No module named loyout


Reading an exiting Layout
~~~~~~~~~~~~~~~~~~~~~~~~~

To read an existing layout it is sufficient to create a Layout object
with, as an argument, a file name with one of the recognized extension.
All files are stored in the ``pstruc['DIRSTRUC']`` directory of the
project. The project root directory is defined in the ``$BASENAME``
environment variable.

.. code:: python

    print pstruc['DIRSTRUC']


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-2-957771f05f1b> in <module>()
    ----> 1 print pstruc['DIRSTRUC']
    

    NameError: name 'pstruc' is not defined


``pstruc`` is a dictionnary which gathers all directories which are used
in ``PyLayers``

.. code:: python

    pstruc


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-3-6d40b18d0400> in <module>()
    ----> 1 pstruc
    

    NameError: name 'pstruc' is not defined


The structure of the ``.osm`` file is shown below

.. code:: python

    %%bash
    
        cd $BASENAME/struc
        ls *.osm


.. parsed-literal::

    DLR.osm
    MADRID-METIS.osm
    MOCAP2.osm
    MOCAP3.osm
    MOCAP-small2.osm
    TA-Office.osm


.. code:: python

    %%bash
        cd $BASENAME/struc
        head DLR.osm
        echo '---'
        tail -17 DLR.osm


.. parsed-literal::

    <?xml version='1.0' encoding='UTF-8'?>
    <osm version='0.6' upload='false' generator='PyLayers'>
    <node id='-212' action='modify' visible='true' lat='47.0100855114' lon='-1.98980710934' />
    <node id='-210' action='modify' visible='true' lat='47.0100789151' lon='-1.9897910381' />
    <node id='-208' action='modify' visible='true' lat='47.0100738861' lon='-1.98977878545' />
    <node id='-206' action='modify' visible='true' lat='47.0100616861' lon='-1.98982814281' />
    <node id='-204' action='modify' visible='true' lat='47.0101583649' lon='-1.98982436917' />
    <node id='-202' action='modify' visible='true' lat='47.0101656174' lon='-1.98981796656' />
    <node id='-200' action='modify' visible='true' lat='47.0101843662' lon='-1.98977935424' />
    <node id='-198' action='modify' visible='true' lat='47.0101791636' lon='-1.98982426816' />
    ---
    <tag k='transition' v='False' />
    </way>
    <way id='-10000123' action='modify' visible='true'>
    <nd ref='-200' />
    <nd ref='-100' />
    <tag k='name' v='WALL' />
    <tag k='z' v="('0.0', '3.0')" />
    <tag k='transition' v='False' />
    </way>
    <way id='-10000124' action='modify' visible='true'>
    <nd ref='-166' />
    <nd ref='-188' />
    <tag k='name' v='WALL' />
    <tag k='z' v="('0.0', '3.0')" />
    <tag k='transition' v='False' />
    </way>
    </osm>


To read a new layout in osm format :

.. code:: python

    L=Layout('DLR.ini')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-6-dc64c84c00f1> in <module>()
    ----> 1 L=Layout('DLR.ini')
    

    NameError: name 'Layout' is not defined


.. code:: python

    fig,ax=L.showGs()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-7-b9ee126fbab0> in <module>()
    ----> 1 fig,ax=L.showGs()
    

    NameError: name 'L' is not defined


.. code:: python

    L.info()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-8-468dc52db9b2> in <module>()
    ----> 1 L.info()
    

    NameError: name 'L' is not defined


The different graphs associated with the layout are then built

.. code:: python

    L.build()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-9-63002b766909> in <module>()
    ----> 1 L.build()
    

    NameError: name 'L' is not defined


The topological graph :math:`\mathcal{G}_t` or graph of non overlapping
cycles.

.. code:: python

    f,a=L.showG('t')
    b=plt.axis('off')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-10-b2e03138f16c> in <module>()
    ----> 1 f,a=L.showG('t')
          2 b=plt.axis('off')


    NameError: name 'L' is not defined


The graph of room :math:`\mathcal{G}_r`. Two rooms which share at least
a wall are connected. Two rooms which share only a corner (punctual
connection) are not connected

.. code:: python

    f,a=L.showG('r')
    b=plt.axis('off')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-11-691b11b3fe3c> in <module>()
    ----> 1 f,a=L.showG('r')
          2 b=plt.axis('off')


    NameError: name 'L' is not defined


The graph of waypath :math:`\mathcal{G}_w`. This graph is used for agent
mobility. This allows to determine the shortest path between 2 rooms.
This information could be included in the osm file. This is not the case
yet

.. code:: python

    f,a=L.showG('w')
    b=plt.axis('off')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-12-17ec2c11acac> in <module>()
    ----> 1 f,a=L.showG('w')
          2 b=plt.axis('off')


    NameError: name 'L' is not defined


The graph of visibility :math:`\mathcal{G_v}`

.. code:: python

    f,a=L.showG('v')
    b=plt.axis('off')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-13-7202a7193478> in <module>()
    ----> 1 f,a=L.showG('v')
          2 b=plt.axis('off')


    NameError: name 'L' is not defined


The graph of interactions :math:`\mathcal{G}_i` used to determine the
ray signatures.

.. code:: python

    f=plt.figure(figsize=(15,15))
    a = f.gca()
    f,a=L.showG('i',fig=f,ax=a)
    b= plt.axis('off')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-14-568749c658a8> in <module>()
    ----> 1 f=plt.figure(figsize=(15,15))
          2 a = f.gca()
          3 f,a=L.showG('i',fig=f,ax=a)
          4 b= plt.axis('off')


    NameError: name 'plt' is not defined


The display options dictionnary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    L.info()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-15-468dc52db9b2> in <module>()
    ----> 1 L.info()
    

    NameError: name 'L' is not defined


The layout can be displayed using matplotlib ploting primitive. Several
display options are specified in the display dictionnary. Those options
are exploited in ``showGs()`` vizualisation method.

.. code:: python

    L.display


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-16-6a8994bd33b4> in <module>()
    ----> 1 L.display
    

    NameError: name 'L' is not defined


Layers
^^^^^^

-  'layer' : list , []
-  'layerset',list, list of available layers
-  'layers', list , []
-  'activelayer', str , 'WINDOW\_GLASS'

-  'alpha', float , 0.5 , overlay transparency
-  'box', tuple , (-20,20,-10,10), (xmin xmax,ymin,ymax)

Strings
^^^^^^^

-  'title' : str , 'Init'
-  'fileoverlay' : str , 'TA-Office.png'

Sizes
^^^^^

-  'fontsize', float , 10
-  'ndsize', float , 10
-  'ndlblsize' : float 20
-  'edlblsize' : float , 20

Booleans
^^^^^^^^

-  'edlabel', boolean, False
-  'ticksoff',boolean, True
-  'scaled' : boolean , True
-  'subseg' : boolean , True
-  'nodes', boolean , True
-  'visu', boolean , False
-  'edges', boolean , True
-  'clear', boolean, False
-  'overlay', boolean , False
-  'thin', boolean , False , If True trace all segments with thickness 1
-  'ndlabel',boolean, If True display node labels
-  'ednodes', boolean, True

Interactive Editor
~~~~~~~~~~~~~~~~~~

The command L.editor() launches an interactive editor. The state machine
is implemented in module ``pylayers.gis.selectl.py``.

To have an idea of all available options, look in the
```pylayers.gis.SelectL`` <http://pylayers.github.io/pylayers/_modules/pylayers/gis/selectl.html#SelectL.new_state>`__
module

All bug correction and ergonomic improvement of this editor is welcome.
Just pull request your modifications.

PyLayers comes along with a low level structure editor based on
``matplotlib`` which can be invoqued using the ``editor()`` method. This
editor is more suited for modyfing constitutive properties of walls. In
the future a dedicated plugin in ``JOSM`` could be a much better
solution.

There are two different modes of edition

-  A create points mode CP

::

    + left clic   : free point
    + right clic  : same x point
    + center clic : same y point

-  A create segments mode

   -  left clic : select point 1
   -  left clic : select point 2
   -  left clic : create a segment between point 1 and point 2

**m** : to switch from one mode to an other

**i** : to return to init state

Image overlay
^^^^^^^^^^^^^

It is useful while editing a layout to have an overlay of an image in
order to help placing points. The image overlay can either be an url or
a filename. In that case the file is stored in

.. code:: python

    L=Layout()
    L.display['fileoverlay']='http://images.wikia.com/theoffice/images/9/9e/Layout.jpg'


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-17-cead63e6ceaa> in <module>()
    ----> 1 L=Layout()
          2 L.display['fileoverlay']='http://images.wikia.com/theoffice/images/9/9e/Layout.jpg'


    NameError: name 'Layout' is not defined


.. code:: python

    L.display['overlay']=True
    L.display['alpha']=1
    L.display['scaled']=False
    L.display['ticksoff']=False
    L.display['inverse']=True


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-18-ee2fbd9090b8> in <module>()
    ----> 1 L.display['overlay']=True
          2 L.display['alpha']=1
          3 L.display['scaled']=False
          4 L.display['ticksoff']=False
          5 L.display['inverse']=True


    NameError: name 'L' is not defined


.. code:: python

    plt.figure(figsize=(10,10))
    L.showGs()


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-19-9bcb9acc34ba> in <module>()
    ----> 1 plt.figure(figsize=(10,10))
          2 L.showGs()


    NameError: name 'plt' is not defined


Scaling the figure overlay
^^^^^^^^^^^^^^^^^^^^^^^^^^

Before going further it is necessary :

-  to place the global origin
-  to precise the vertical and horizontal scale of the image

This is done by the following commands :

-  'i' : back to init state
-  'm' : goes to CP state
-  'o' : define the origin
-  'left click' on the point of the figure chasen as the origin
-  'left click' on a point at a known distance from the origin along x
   axis. Fill the dialog box with the actual distance (expressed in
   meters) between the two points.
-  'left click' on a point at a known distance from the origin along y
   axis. Fill the dialog box with the actual distance (expressed in
   meters) between the two points.

In that sequence of operation it is useful to rescale the figure with
'r'.

At that stage, it is possible to start creating points

::

        'b'  : selct a segment
        'l'  : select activelayer
        'i'  : back to init state
        'e'  : edit segment
        't'  : translate  structure
        'h'  : add subsegment
        'd'  : delete subsegment
        'r'  : refresh
        'o'  : toggle overlay
        'm'  : toggle mode (point or segment)
        'z'  : change display parameters
        'q'  : quit interactive mode
        'x'  : save .str2 file
        'w'  : display all layers

Vizualisation of the layout
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    L = Layout('TA-Office.ini')
    L.dumpr()
    fig = plt.figure(figsize=(25,25))
    ax = fig.gca()
    fig,ax = L.showG(fig=fig,ax=ax,graph='s',labels=True,font_size=9,node_size=220,node_color='c')
    a = plt.axis('off')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-20-19a72aa609b6> in <module>()
    ----> 1 L = Layout('TA-Office.ini')
          2 L.dumpr()
          3 fig = plt.figure(figsize=(25,25))
          4 ax = fig.gca()
          5 fig,ax = L.showG(fig=fig,ax=ax,graph='s',labels=True,font_size=9,node_size=220,node_color='c')


    NameError: name 'Layout' is not defined


Each node of :math:`\mathcal{G}_s` with a negative index is a point.

Each node of :math:`\mathcal{G}_s` with a positive index corresponds to
a segment (wall,door,window,...).

The segment name is the key of the **slab** dictionnary.

`Multi Subsegments <./Multisubsegments.ipynb>`__
