
.. code-block:: python

    from pylayers.gis.layout import *
    from IPython.display import Image
    import os


.. parsed-literal::

    <matplotlib.figure.Figure at 0x41f3f10>


The ``ls()`` method
-------------------


The ``ls()`` method lists the layout file which are available in the
``struc`` directory of the current project.

.. code-block:: python

    L=Layout()
    L.ls()



.. parsed-literal::

    ['11D-E1.ini',
     'DLR.ini',
     'DLR2.ini',
     'Lstruc.ini',
     'TA-Office.ini',
     'W2PTIN.ini',
     'WHERE1.ini',
     'defstr.ini',
     'defstr3.ini',
     'klepal.ini',
     'w2ptin.ini',
     'where1.ini']



.. code-block:: python

    L.ls(typ='osm')



.. parsed-literal::

    ['DLR.osm', 'DLR2.osm', 'where1.osm']



.. code-block:: python

    L=Layout('DLR.osm')

.. code-block:: python

    L.Gs.node[1]



.. parsed-literal::

    {'name': 'WALL',
     'norm': array([-0.52106509, -0.85351694,  0.        ]),
     'transition': False,
     'z': ('0.0', '3.0')}



.. code-block:: python

    L.saveosm('DLR.osm')
.. code-block:: python

    Image('../../data/struc/images/TA-Office.png')



.. image:: DEV-Layout_files/DEV-Layout_8_0.png



.. code-block:: python

    fig,ax=L.showGs()


.. image:: DEV-Layout_files/DEV-Layout_9_0.png


.. code-block:: python

    L=Layout('TA-Office.ini')
.. code-block:: python

    L.showGs()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x42ff6d0>,
     <matplotlib.axes.AxesSubplot at 0x3484fd0>)




.. image:: DEV-Layout_files/DEV-Layout_11_1.png


.. code-block:: python

    L=Layout('11D-E1.ini')
    L



.. parsed-literal::

    
    ----------------
    11D-E1.ini
    ----------------
    
    Number of points  : 214
    Number of segments  : 199
    Number of sub segments  : 0
    Number of cycles  : 0
    Number of rooms  : 0
    degree 0 : []
    degree 1 : [  -1 -214 -213 -212 -211 -210 -206 -205 -204 -202 -201 -200 -198 -197 -196
     -195 -194 -193 -191 -188 -184 -183 -182 -180 -179 -178 -177 -175 -174 -172
     -171 -167 -166 -164 -162 -161 -159 -156 -154 -152 -150 -148 -147 -145 -144
     -142 -140 -139 -137 -135 -134 -132 -130 -129 -128 -126 -124 -123 -122 -120
     -118 -116 -115 -113 -112 -110 -109 -107 -105 -103 -102 -100  -98  -96  -94
      -92  -90  -88  -87  -85  -83  -81  -80  -78  -76  -75  -73  -71  -70  -68
      -67  -65  -63  -62  -59  -57  -55  -54  -52  -51  -49  -48  -46  -43  -41
      -40  -39  -37  -36  -33  -31  -30  -28  -26  -25  -24  -21  -20  -18  -17
      -15  -13  -12  -10   -9   -8   -7   -6   -5   -4   -3   -2]
    degree 2 : 9
    degree 3 : 44
    degree 4 : [-190 -185 -181 -176 -173 -170 -163 -160 -149 -143 -138 -131 -127 -121 -114
     -108 -101  -82  -79  -72  -66  -53  -50  -38  -35  -22  -19  -14  -11]
    
    xrange :(-10.462, 73.369)
    yrange :(-0.096, 15.004)
    
    Useful dictionnaries
    ----------------
    sl {slab name : slab dictionary}
    name :  {slab :seglist} 
    
    Useful arrays
    ----------------
    tsg : get segment index in Gs from tahe
    tgs : get segment index in tahe from Gs
    lsss : list of segments with sub-segment
    sla : associated slab name
    stridess : stride for adressing sub segment 
    degree : degree of nodes 




.. code-block:: python

    L.showG('s',figsize=(20,10))



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x293b850>,
     <matplotlib.axes.AxesSubplot at 0x61d6110>)




.. image:: DEV-Layout_files/DEV-Layout_13_1.png


.. code-block:: python

    L=Layout('klepal.ini')
    L



.. parsed-literal::

    
    ----------------
    klepal.ini
    Image('/home/uguen/Bureau/P1/struc/images/IMG-Layout-Klepal.png')
    ----------------
    
    Number of points  : 25
    Number of segments  : 11
    Number of sub segments  : 0
    Number of cycles  : 0
    Number of rooms  : 0
    degree 0 : [-25 -24 -23 -22 -21 -20 -19 -13  -5  -4  -3  -1]
    degree 1 : [ -2 -18 -17 -16 -15 -14]
    degree 2 : 5
    degree 3 : 2
    
    xrange :(-17.872, 18.155)
    yrange :(-6.869, 6.51)
    
    Useful dictionnaries
    ----------------
    sl {slab name : slab dictionary}
    name :  {slab :seglist} 
    
    Useful arrays
    ----------------
    tsg : get segment index in Gs from tahe
    tgs : get segment index in tahe from Gs
    lsss : list of segments with sub-segment
    sla : associated slab name
    stridess : stride for adressing sub segment 
    degree : degree of nodes 




This Layout is still in construction

.. code-block:: python

    L.showGs()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f6271f91b10>,
     <matplotlib.axes.AxesSubplot at 0x7f6271f905d0>)




.. image:: DEV-Layout_files/DEV-Layout_16_1.png


.. code-block:: python

    L=Layout('W2PTIN.ini')
.. code-block:: python

    L.showGs()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f6271f59f10>,
     <matplotlib.axes.AxesSubplot at 0x61bebd0>)




.. image:: DEV-Layout_files/DEV-Layout_18_1.png


.. code-block:: python

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


