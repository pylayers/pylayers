How to edit a new structure ? 
=============================

Terminology
-----------

Terms **nodes** and **edges** are used when refering to graph elements of :math:`\mathcal{G}_{x}(\mathcal{V},\mathcal{E})`. 
Terms **points** and **segments** are used when refering to geometrical object. 

The `.str2` file 
----------------


An indoor structure is originally described in a very simple text format. 
This ASCII file format is the simplest description of the layout. 

A `str2` ASCII file is a list of points with their coordinates and a list of segments associating 
those points by pairs. 

The description of the layout exploits a **non overlaping  rule**,

(which is ( `or should be` ) checked regularly when creating a new layout for not braking the consistency of the associated graph description).

.. warning::

    **Non overlaping rule**:
    A segment must never overlap any point of the layout  

This rule aims maintaining, consistent spatial relationships between points and segments. 
For example if there is a long wall with many doors, each of the short successive segments
must be explicitely described. It is not authorized to model long walls by
a long segment which overlaps all the doors. 



Creating the layout
-------------------

.. code:: python 
    
    from pylayers.gis.layout import *

    L = Layout('fname.str2')

Either the file `fname.str2` exists in which case it is loaded, otherwise a
void structure is created. 


Layout editor matplotlib GUI 
----------------------------

It exist a very simple structure editor based on matplotlib events handling. 
This editor is simply as a particular method of the Layout object

.. code:: python 

     L.editor()

At that state a matplolib figure is opened which can interact either with
mouse or keayboard events. 


Points creation mode  
--------------------

The point creation mode is reachable through the key `m`
Once in point edition mode, the clouse click is used to create a new point in
the layout.

+ left click : add a free point 
+ right click : add a new point with same x coordinate as the previously added
  point 
+ center click : add a new point with same y coordinate as the previous added
  point 

For the point edition phase it is generally convenient to have an image
overlay. This is possible using the `o` key which toggle the current overlay
image.  

To delete a single point, select the point and type `d`. 
To delete a group of point use the zoom to delimitate the region to delete and
type the `Suppr` key. 


Walls edition mode
------------------


Modification of display parameters
----------------------------------

Typing key `z` open a window box for modifying various parameters

  - filename 
  - fileoverlay 
  - box 
  - alpha 


List of keyboard shortcuts 
--------------------------

+ `x` : save the current layout 
+ `i` : back in init mode  
+ `q` : quit intercative mode 
+ `m` : toggle point creation mode
+ `z` : modify display parameters  



Using graphs
------------

The key design concept of the ray tracing tool is to pre-process
the layout description exploiting the **graph** abstract data structure
exploiting the `networkx <http:// >` package. z


see :ref:`Uguen2012`_


