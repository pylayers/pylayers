
.. code:: python

    from pylayers.gis.ezone import *
    %matplotlib inline
An Ezone is an objet which gathers vector and raster information from a
1 degree by 1 degree tile of earth

.. code:: python

    prefix=enctile(-1.5,10.5)
    print prefix

.. parsed-literal::

    N10W002


.. code:: python

    prefix



.. parsed-literal::

    'N10W002'



.. code:: python

    dectile(prefix=prefix)



.. parsed-literal::

    (-2, -1, 10, 11)



.. code:: python

    int('08')



.. parsed-literal::

    8



.. code:: python

    E=Ezone(prefix)
.. code:: python

    E.prefix



.. parsed-literal::

    'N10W002'



An ``Ezone`` can be obtained from a point (longitude,Latitude)

.. code:: python

    r=E.getdem()

.. parsed-literal::

    Download srtm file
    SRTMDownloader - server= dds.cr.usgs.gov, directory=srtm/version2_1/SRTM3.
    no aster file for this point


.. code:: python

    E.saveh5()
.. code:: python

    f,a = E.show(source='srtm',clim=[0,500])


.. image:: Downloading_files/Downloading_11_0.png


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


