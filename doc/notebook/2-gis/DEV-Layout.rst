
.. code-block:: python

    from pylayers.gis.layout import *
    from IPython.display import Image
    import os


.. parsed-literal::

    <matplotlib.figure.Figure at 0xa34cf4c>


The ``ls()`` method
-------------------

The ``ls()`` method lists the layout file which are available in the
``struc`` directory of the current project.

.. code-block:: python

    L=Layout()
    L.ls()



.. parsed-literal::

    ['DLR.ini',
     'DLR2.ini',
     'TA-Office.ini',
     'W2PTIN.ini',
     'WHERE1.ini',
     'WHERE1_old.ini',
     'defstr3.ini']



.. code-block:: python

    L.ls(typ='osm')



.. parsed-literal::

    ['DLR.osm']



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

    (<matplotlib.figure.Figure at 0xac84f6c>,
     <matplotlib.axes.AxesSubplot at 0xac8404c>)




.. image:: DEV-Layout_files/DEV-Layout_11_1.png


.. code-block:: python

    L=Layout('11D-E1.ini')
    L



.. parsed-literal::

    <repr(<pylayers.gis.layout.Layout at 0xa8e398c>) failed: AttributeError: 'Layout' object has no attribute 'degree'>



.. code-block:: python

    L.showG('s',figsize=(20,10))



.. parsed-literal::

    (<matplotlib.figure.Figure at 0xa8c052c>,
     <matplotlib.axes.AxesSubplot at 0xac84aac>)




.. image:: DEV-Layout_files/DEV-Layout_13_1.png


.. code-block:: python

    L=Layout('klepal.ini')
    L



.. parsed-literal::

    <repr(<pylayers.gis.layout.Layout at 0xa47354c>) failed: AttributeError: 'Layout' object has no attribute 'degree'>



This Layout is still in construction

.. code-block:: python

    L.showGs()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0xa53d78c>,
     <matplotlib.axes.AxesSubplot at 0xa543c0c>)




.. image:: DEV-Layout_files/DEV-Layout_16_1.png


.. code-block:: python

    L=Layout('W2PTIN.ini')
.. code-block:: python

    L.showGs()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0xa7cff6c>,
     <matplotlib.axes.AxesSubplot at 0xa7cf6ec>)




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


