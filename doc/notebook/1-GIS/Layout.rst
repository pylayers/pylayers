
Description of the propagation environment
==========================================

The ``Layout`` class contains the data structure for describing an
Indoor environment. It implements the different graphs helping the
implementation of the ray tracing. The class is implemented in the
```layout.py`` <http://pylayers.github.io/pylayers/modules/pylayers.gis.layout.html>`__
module.

.. code:: python

    from pylayers.gis.layout import *
    from IPython.display import Image
    import os
    %matplotlib inline


.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


Getting the list of all available Layouts : the ``ls()`` method
---------------------------------------------------------------

Creating a default Layout is as simple as :

.. code:: python

    L=Layout()
    L




.. parsed-literal::

    
    ----------------
    defstr.ini
    ----------------
    
    Number of points  : 8
    Number of segments  : 9
    Number of sub segments  : 3
    Number of cycles  : 3
    Number of rooms  : 2
    degree 0 : []
    degree 1 : [-8 -7]
    number of node point of degree 2 : 4
    number of node point of degree 3 : 2
    
    xrange :(758.49, 768.516)
    yrange :(1111.9, 1115.963)
    
    Useful dictionnaries
    ----------------
    dca {cycle : []} cycle with an airwall
    sl {slab name : slab dictionary}
    name :  {slab :seglist} 
    
    Useful arrays
    ----------------
    pt : numpy array of points 
    normal : numpy array of normal 
    offset : numpy array of offset 
    tsg : get segment index in Gs from tahe
    isss :  sub-segment index above Nsmax
    tgs : get segment index in tahe from Gs
    lsss : list of segments with sub-segment
    sla : list of all slab names (Nsmax+Nss+1)
    degree : degree of nodes 



Querying the default file name as simple as :

.. code:: python

    L.filename




.. parsed-literal::

    'defstr.ini'



The ``ls()`` method lists the layout files which are available in the
``struc`` directory of your current project, which is set up via the
$BASENAME environment variable which is crucial to be early defined in
order PyLayers find its way to the good directories. Over the
development process, the layout data format has evolved quite a lot, the
most simple is an ``ini`` key-value text file.

.. code:: python

    L.ls('ini')




.. parsed-literal::

    ['CORM1.ini',
     'DLR.ini',
     'DLR2.ini',
     'MADRID-METIS.ini',
     'MOCAP-small.ini',
     'MOCAP-small2.ini',
     'MOCAP-small3.ini',
     'MOCAP.ini',
     'MOCAPext.ini',
     'Scene.ini',
     'TA-Office.ini',
     'TA-OfficeAir.ini',
     'W2PTIN.ini',
     'WHERE1.ini',
     'WHERE2.ini',
     'd24.ini',
     'defstr.ini',
     'defstr3.ini',
     'homeK_vf.ini',
     'klepal.ini',
     'nicta.ini',
     'scat1.ini',
     'scat2.ini',
     'scattering.ini',
     'test.ini']



.. code:: python

    L=Layout('DLR.ini')

.. code:: python

    f,a=L.showG('s')



.. image:: Layout_files/Layout_11_0.png


To check which are the used slabs :

.. code:: python

    Slabs = np.unique(L.sla)
    for s in Slabs:
        if s in L.sl:
            print L.sl[s]


.. parsed-literal::

    3D_WINDOW_GLASS : GLASS | AIR | GLASS | [0.005, 0.005, 0.005]
    
    AIR : AIR | [0.02]
    
    DOOR : WOOD | [0.03]
    
    METAL : METAL | [0.1]
    
    PARTITION : PLASTER | [0.1]
    
    WALL : BRICK | [0.07]
    


Let's load an other layout

.. code:: python

    L=Layout('WHERE1.ini')

The showG method provides many vizualization of the layout

.. code:: python

    f,a=L.showG('s',airwalls=False,figsize=(20,10))



.. image:: Layout_files/Layout_17_0.png


.. code:: python

    L=Layout('W2PTIN.ini')


::


    ---------------------------------------------------------------------------

    AssertionError                            Traceback (most recent call last)

    <ipython-input-10-366aeaf5fde2> in <module>()
    ----> 1 L=Layout('W2PTIN.ini')
    

    /home/uguen/Documents/rch/devel/pylayers/pylayers/gis/layout.pyc in __init__(self, _filename, _filematini, _fileslabini, _filefur, force, check)
        411         # check layout integrity (default)
        412         if check:
    --> 413             self.check()
        414         self.boundary()
        415 


    /home/uguen/Documents/rch/devel/pylayers/pylayers/gis/layout.pyc in check(self, level)
        719             deg0 = filter(lambda x: nx.degree(self.Gs,x)==0,upnt)
        720             deg1 = filter(lambda x: nx.degree(self.Gs,x)==1,upnt)
    --> 721             assert (len(deg0)==0), "It exists degree 0 points :  %r" % deg0
        722             assert (len(deg1)==0), "It exists degree 1 points : %r" % deg1
        723 


    AssertionError: It exists degree 0 points :  [-110, -109, -108, -103]


.. code:: python

    f,a = L.showG('s')



.. image:: Layout_files/Layout_19_0.png


The useful numpy arrays of the Layout
-------------------------------------

The layout data structure is a mix between graph and numpy array. numpy
arrays are used when high performance is required while graph structure
is convenient when dealing with different specific tasks. The tricky
thing for the mind is to have to transcode between node index excluding
0 and numpy array index including 0. Below are listed various useful
numpy array which are mostly used internally.

-  tsg : get segment index in Gs from tahe
-  isss : sub-segment index above Nsmax
-  tgs : get segment index in tahe from Gs
-  lsss : list of segments with sub-segment
-  sla : list of all slab names (Nsmax+Nss+1)
-  degree : degree of nodes

``pt`` the array of points
~~~~~~~~~~~~~~~~~~~~~~~~~~

The point coordinates are stored in two different places (which in
principle is a bad thing to do !).

::

    L.Gs.pos : in a dictionnary form (key is the point negative index)
    L.pt : in a numpy array

.. code:: python

    print np.shape(L.pt)
    print len(filter(lambda x: x<0,L.Gs.pos))


.. parsed-literal::

    (2, 278)
    278


This dual storage is chosen (temporarily ? ) for computational
efficiency reason. The priority goes to the graph and the numpy array is
calculated at the end of the edition in the ``Layout.g2npy`` method
(graph to numpy) which is in charge of the conversion.

tahe (tail-head)
~~~~~~~~~~~~~~~~

``tahe`` is a :math:`(2\times N_{s})` where :math:`N_s` denotes the
number of segment. The first line is the tail index of the segment
:math:`k` and the second line is the head of the segment :math:`k`.
Where :math:`k` is the index of a given segment (starting in 0).

.. code:: python

    L.build()

The figure below illustrates a Layout and a surimposition of the graph
of cycles :math:`\mathcal{G}_c`. Those cycles are automatically
extracted from a well defined layout. This concept of **cycles** is
central in the ray determination algorithm which is implemented in
PyLayers. Notice that the exterior region is the cycle indexed by 0. All
the rooms which have a common frontier with the exterior cycle are here
connected to the origin (corresponding to exterior cycle).

.. code:: python

    f,a = L.showG('s')
    nx.draw(L.Gc,L.Gc.pos)



.. image:: Layout_files/Layout_32_0.png


.. code:: python

    nx.draw_networkx_nodes(L.Gi,L.Gi.pos,node_color='blue',node_size=1)
    nx.draw_networkx_edges(L.Gi,L.Gi.pos,node_color='blue',node_size=1)




.. parsed-literal::

    <matplotlib.collections.LineCollection at 0x2b1d75935950>




.. image:: Layout_files/Layout_33_1.png


``tgs`` : trancodage from graph indexing to numpy array indexing
----------------------------------------------------------------

``tgs`` is an array with length :math:`N_s`\ +1. The index 0 is not used
because none segment has 0 as an index.

.. code:: python

    ns = 5
    utahe = L.tgs[ns]

.. code:: python

    tahe =  L.tahe[:,utahe]

.. code:: python

    ptail = L.pt[:,tahe[0]]
    phead = L.pt[:,tahe[1]]

.. code:: python

    print ptail


.. parsed-literal::

    [-28.081  10.923]


.. code:: python

    print phead


.. parsed-literal::

    [-28.118  14.857]


.. code:: python

    L.Gs.node[5]




.. parsed-literal::

    {'connect': [-286, -292],
     'name': 'CONCRETE_20CM3D',
     'ncycles': [6, 0],
     'norm': array([-0.99995577, -0.00940477,  0.        ]),
     'offset': 0,
     'ss_name': ['3D_WINDOW_GLASS'],
     'ss_offset': [0],
     'ss_z': [(1.5, 2.5)],
     'transition': False,
     'z': (0.0, 3.0)}



.. code:: python

    print L.Gs.pos[-8]
    print L.Gs.pos[-139]


.. parsed-literal::

    (31.687, 11.252)
    (5.037, 10.963)


.. code:: python

    aseg = np.array([4,7,134])

.. code:: python

    print np.shape(aseg)


.. parsed-literal::

    (3,)


.. code:: python

    pt  = L.tahe[:,L.tgs[aseg]][0,:]
    ph = L.tahe[:,L.tgs[aseg]][1,:]
    pth = np.vstack((pt,ph))

.. code:: python

    np.shape(pth)




.. parsed-literal::

    (2, 3)



``Layout.seg2pts`` a function for getting points coordinates from segment number array
--------------------------------------------------------------------------------------

.. code:: python

    L.seg2pts(aseg)




.. parsed-literal::

    array([[-28.081, -27.833,   0.454],
           [ 10.923,  10.686,   4.805],
           [-27.836, -27.835,   0.457],
           [ 10.926,  10.891,   4.529]])



.. code:: python

    aseg = array(filter(lambda x: x>0,L.Gs.nodes()))
    pth = L.seg2pts(aseg)

.. code:: python

    from pylayers.util.plotutil import displot

.. code:: python

    displot(pth[0:2,:],pth[2:,:])
    plt.axis('off')




.. parsed-literal::

    (-30.0, 40.0, 4.0, 18.0)




.. image:: Layout_files/Layout_51_1.png

