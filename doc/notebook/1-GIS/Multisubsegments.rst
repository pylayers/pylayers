
Effect of Modyfiying the Nature of Sub-Segments
===============================================

.. code:: python

    from pylayers.simul.link import *
    from pylayers.antprop.rays import *
    from pylayers.gis.layout import *
    from pylayers.antprop.signature import *
    import pylayers.signal.bsignal as bs
    import pylayers.signal.waveform as wvf 
    from pylayers.simul.simulem import *
    import matplotlib.pyplot as plt 
    %matplotlib inline

.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


This section presents a simple Ray Tracing simulation with different
material properties of a subsegment separating 2 rooms.

.. code:: python

    L=Layout('defstr3.ini')
    f,a=L.showG('s',subseg=True,figsize=(10,10))


.. image:: Multisubsegments_files/Multisubsegments_3_0.png


The studied configuration is composed of a simple 2 rooms building
separated by a subsegment which has a multi subsegment attribute. The
attribute of the subsegment can be changed with the method
```chgmss`` <http://pylayers.github.io/pylayers/modules/generated/pylayers.gis.layout.Layout.chgmss.html>`__
(change multisubsegment)

.. code:: python

    L.chgmss(1,ss_name=['WOOD','AIR','WOOD'],ss_z =[(0.0,2.7),(2.7,2.8),(2.8,3)],ss_offset=[0,0,0])
As the Layout structure has been modified it is required to rebuild the
structure.

.. code:: python

    L.build()
    L.save()

.. parsed-literal::

    structure saved in  defstr3.str2
    structure saved in  defstr3.ini


The :math:`\mathcal{G}_s` graph dictionnary has the following structure

.. code:: python

    L.Gs.node



.. parsed-literal::

    {-8: {},
     -7: {},
     -6: {},
     -5: {},
     -4: {},
     -3: {},
     -2: {},
     -1: {},
     1: {'connect': [-8, -7],
      'name': 'PARTITION',
      'ncycles': [1, 2],
      'norm': array([-0.999982  , -0.00599989,  0.        ]),
      'offset': 0,
      'ss_name': ['WOOD', 'AIR', 'WOOD'],
      'ss_offset': [0, 0, 0],
      'ss_z': [(0.0, 2.7), (2.7, 2.8), (2.8, 3)],
      'transition': True,
      'z': (0.0, 3.0)},
     2: {'connect': [-8, -2],
      'name': 'WALL',
      'ncycles': [1, 2],
      'norm': array([ 0.99997778,  0.00666652,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     3: {'connect': [-7, -5],
      'name': 'WALL',
      'ncycles': [1, 2],
      'norm': array([-0.99997775, -0.00667097,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     4: {'connect': [-6, -1],
      'name': 'WALL',
      'ncycles': [2, 0],
      'norm': array([ 0.99997888,  0.00649986,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     5: {'connect': [-6, -5],
      'name': 'WALL',
      'ncycles': [2, 0],
      'norm': array([-0.00619988,  0.99998078,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     6: {'connect': [-5, -4],
      'name': 'WALL',
      'ncycles': [1, 0],
      'norm': array([-0.00639987,  0.99997952,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     7: {'connect': [-4, -3],
      'name': 'WALL',
      'ncycles': [1, 0],
      'norm': array([ 0.99997887,  0.00650149,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     8: {'connect': [-3, -2],
      'name': 'WALL',
      'ncycles': [1, 0],
      'norm': array([ 0.00639987, -0.99997952,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)},
     9: {'connect': [-2, -1],
      'name': 'WALL',
      'ncycles': [2, 0],
      'norm': array([ 0.00639987, -0.99997952,  0.        ]),
      'offset': 0,
      'transition': False,
      'z': (0.0, 3.0)}}





We define now two points which are the termination of a radio link.

.. code:: python

    tx=np.array([759,1114,1.0])
    rx=np.array([767,1114,1.5])
.. code:: python

    L.chgmss(1,ss_name=['WOOD','AIR','WOOD'],ss_z =[(0.0,2.7),(2.7,2.8),(2.8,3)],ss_offset=[0,0,0])
    L.save()
    Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))

.. parsed-literal::

    structure saved in  defstr3.str2
    structure saved in  defstr3.ini


.. code:: python

    Lk



.. parsed-literal::

    filename: Links_0_defstr3.ini.h5
    Link Parameters :
    ------- --------
    Layout : defstr3.ini
    
    Node a   
    ------  
    position : [  7.59000000e+02   1.11400000e+03   1.00000000e+00]
    Antenna : Omni
    Rotation matrice : 
     [[ 1.  0.  0.]
     [ 0.  1.  0.]
     [ 0.  0.  1.]]
    
    Node b   
    ------  
    position : [  767.   1114.      1.5]
    Antenna : Omni
    Rotation matrice : 
     [[ 1.  0.  0.]
     [ 0.  1.  0.]
     [ 0.  0.  1.]]
    
    Link evaluation information : 
    ----------------------------- 
    distance :  8.016 m 
    delay : 26.719 ns
    fmin (fGHz) : 2.0
    fmax (fGHz) : 11.0
    fstep (fGHz) : 0.05
     



.. code:: python

    f,a=Lk.show()


.. image:: Multisubsegments_files/Multisubsegments_15_0.png


On the figure above, we can see the Tx and Rx each placed in a different
room appart from a wall with a subsegement placed in the middle. Then
for evaluating the radio link, simply type:

.. code:: python

    ak,tauk=Lk.eval(force=True)

.. parsed-literal::

    Signatures'> from 2_1_3 saved
    Rays'> from 3_0_1 saved
    Ctilde'> from 0_1_0 saved
    Tchannel'> from 0_1_0_0_0_1_1 saved


At that point the channel has been evaluated and all the data stored in
an ``hdf5`` file

The different members of the link are
-------------------------------------

The Signature of the radio channel is in ``Lk.Si``, the 3D rays are in
``Lk.R``, the propagation channel is in ``Lk.C`` and the transmission
channel is in ``Lk.H``

.. code:: python

    Lk.R



.. parsed-literal::

    Rays3D
    ----------
    1 / 1 : [0]
    2 / 6 : [1 2 3 4 5 6]
    3 / 18 : [ 7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]
    4 / 37 : [25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
     50 51 52 53 54 55 56 57 58 59 60 61]
    5 / 42 : [ 62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
      80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97
      98  99 100 101 102 103]
    6 / 26 : [104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121
     122 123 124 125 126 127 128 129]
    -----
    ni : 581
    nl : 1292




.. code:: python

    Lk.C



.. parsed-literal::

    Ctilde
    ---------
    (130, 181)
    Nray : 130
    fmin(GHz) : 2.0
    fmax(GHz): 11.0
    Nfreq : 181



.. code:: python

    f = plt.figure(figsize=(10,10))
    f,a=Lk.C.show(cmap='jet',fig=f,typ='l20')


.. image:: Multisubsegments_files/Multisubsegments_23_0.png


.. code:: python

    Lk.H



.. parsed-literal::

    freq :2.0 10.0 161
    shape  :(130, 161)
    tau :26.7186992365 73.2944728109
    dist :8.01560977094 21.9883418433
    Friis factor -j c/(4 pi f) has been applied



.. code:: python

    fGHz=np.arange(2,6,0.5)
    wav = wvf.Waveform(fcGHz=4,bandGHz=1.5)
    wav.show()


.. image:: Multisubsegments_files/Multisubsegments_25_0.png


.. code:: python

    cir = Lk.H.applywavB(wav.sf)
.. code:: python

    cir.plot(typ=['v'],xmin=20,xmax=80)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f392787d950>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f39275ebfd0>]], dtype=object))




.. image:: Multisubsegments_files/Multisubsegments_27_1.png


.. code:: python

    layer = ['AIR','AIR','AIR']
    L.chgmss(1,ss_name=layer)
    L.Gs.node[1]['ss_name']=layer
    L.g2npy()
    L.save()
    Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))
    Lk.eval(force=True)
    cirair = Lk.H.applywavB(wav.sf)
    cirair.plot(typ=['v'],xmin=20,xmax=80)

.. parsed-literal::

    structure saved in  defstr3.str2
    structure saved in  defstr3.ini
    Signatures'> from 2_1_3 saved
    Rays'> from 3_0_1 saved
    Ctilde'> from 0_1_0 saved
    Tchannel'> from 0_1_0_0_0_1_1 saved




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f39272475d0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f3926fbb690>]], dtype=object))




.. image:: Multisubsegments_files/Multisubsegments_28_2.png


.. code:: python

    layer = ['PARTITION','PARTITION','PARTITION']
    L.chgmss(1,ss_name=layer)
    L.Gs.node[1]['ss_name']=layer
    L.g2npy()
    L.save()
    Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))
    Lk.eval(force=True)
    cirpart = Lk.H.applywavB(wav.sf)
    cirpart.plot(typ=['v'],xmin=20,xmax=80)

.. parsed-literal::

    structure saved in  defstr3.str2
    structure saved in  defstr3.ini
    Signatures'> from 2_1_3 saved
    Rays'> from 3_0_1 saved
    Ctilde'> from 0_1_0 saved
    Tchannel'> from 0_1_0_0_0_1_1 saved




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f392b898d90>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f392b8988d0>]], dtype=object))




.. parsed-literal::

    <matplotlib.figure.Figure at 0x7f392b1aa2d0>



.. image:: Multisubsegments_files/Multisubsegments_29_3.png


.. code:: python

    layer = ['METAL','METAL','METAL']
    L.chgmss(1,ss_name=layer)
    L.Gs.node[1]['ss_name']=layer
    L.g2npy()
    L.save()
    Lk = DLink(L=L,a=tx,b=rx,Aa=Antenna('Omni'),Ab=Antenna('Omni'))
    Lk.eval(force=True)
    cirmet = Lk.H.applywavB(wav.sf)
    cirmet.plot(typ=['v'],xmin=20,xmax=80)

.. parsed-literal::

    structure saved in  defstr3.str2
    structure saved in  defstr3.ini
    Signatures'> from 2_1_3 saved
    Rays'> from 3_0_1 saved
    Ctilde'> from 0_1_0 saved
    Tchannel'> from 0_1_0_0_0_1_1 saved




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f392793d690>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f3926fb4d10>]], dtype=object))




.. image:: Multisubsegments_files/Multisubsegments_30_2.png


.. code:: python

    #fig2=plt.figure()
    f,a=cirair.plot(typ=['l20'],color='b')
    plt.axis([0,120,-120,-40])
    plt.title('A simple illustration of shadowing effect')
    plt.legend(['air'])
    f,a=cirpart.plot(typ=['l20'],color='k')
    plt.axis([0,120,-120,-40])
    plt.legend(['wood'])
    f,a=cirmet.plot(typ=['l20'],color='r')
    plt.axis([0,120,-120,-40])
    plt.legend(['metal'])



.. parsed-literal::

    <matplotlib.legend.Legend at 0x7f3926d68bd0>




.. image:: Multisubsegments_files/Multisubsegments_31_1.png



.. image:: Multisubsegments_files/Multisubsegments_31_2.png



.. image:: Multisubsegments_files/Multisubsegments_31_3.png


| We have modified successively the nature of the 3 surfaces in the sub
segment placed in the sepataion partition. The first was AIR, the second
WOOD and the third METAL. As the subsegment is placed on the LOS path
the blockage effect is clearly visible.
| The chosen antennas were omnidirectional ``Antenna('Omni')``

.. code:: python

    from IPython.core.display import HTML
    
    def css_styling():
        styles = open("../styles/custom.css", "r").read()
        return HTML(styles)
    css_styling()



.. raw:: html

    <style>
        @font-face {
            font-family: "Computer Modern";
            src: url('http://mirrors.ctan.org/fonts/cm-unicode/fonts/otf/cmunss.otf');
        }
        div.cell{
            width:800px;
            margin-left:16% !important;
            margin-right:auto;
        }
        h1 {
            font-family: Helvetica, serif;
        }
        h4{
            margin-top:12px;
            margin-bottom: 3px;
           }
        div.text_cell_render{
            font-family: Computer Modern, "Helvetica Neue", Arial, Helvetica, Geneva, sans-serif;
            line-height: 145%;
            font-size: 130%;
            width:800px;
            margin-left:auto;
            margin-right:auto;
        }
        .CodeMirror{
                font-family: "Source Code Pro", source-code-pro,Consolas, monospace;
        }
        .prompt{
            display: None;
        }
        .text_cell_render h5 {
            font-weight: 300;
            font-size: 22pt;
            color: #4057A1;
            font-style: italic;
            margin-bottom: .5em;
            margin-top: 0.5em;
            display: block;
        }
        
        .warning{
            color: rgb( 240, 20, 20 )
            }  
    </style>
    <script>
        MathJax.Hub.Config({
                            TeX: {
                               extensions: ["AMSmath.js"]
                               },
                    tex2jax: {
                        inlineMath: [ ['$','$'], ["\\(","\\)"] ],
                        displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
                    },
                    displayAlign: 'center', // Change this to 'center' to center equations.
                    "HTML-CSS": {
                        styles: {'.MathJax_Display': {"margin": 4}}
                    }
            });
    </script>



