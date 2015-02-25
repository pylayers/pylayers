
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



Desciption of antennas
======================

PyLayers has a very rich set of tools for handling antenna radiation
pattern. Antennas can be described in different manners and read from
different specific file formats.

The description goes from a simple antenna gain formula to a full
polarization description ,compressed or not, using scalar or vector
spherical harmonics decomposition.

In the following some features of the ``Antenna`` class are illustrated.
The ``Antenna`` class is stored in the
`antenna.py <http://pylayers.github.io/pylayers/modules/pylayers.antprop.antenna.html>`__
module which is placed in the ``antprop`` module.

.. code:: python

    from pylayers.antprop.antenna import *
    %matplotlib inline

.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


An antenna object can be loaded in specifying an existing antenna file
name as argument of the constructor. Lets start by loading an antenna
from a ``vsh3`` file which correspond to a vector spherical harmonics
representation of an antenna measured in SATIMO near field chamber.

.. code:: python

    A = Antenna('S1R1.vsh3')
    A



.. parsed-literal::

    FileName : S1R1.vsh3
    -----------------------
    fmin : 0.80GHz
    fmax : 5.95GHz
    step : 50.00MHz
    Nf : 104
    Not evaluated



At loading time the antenna is not evaluated. It means that there is not
internally any instanciation of the pattern for a set of angular and
frequency values.

To list all the availble antenna files in the dedicated directory of the
project it is possible to invoke the ``ls()`` method.

Antenna files are supposed to be stored in the sub directory ``ant`` of
the current project. The current project is located with the $BASENAME
environment variable.

.. code:: python

    !echo $BASENAME

.. parsed-literal::

    /home/uguen/Bureau/P1


.. code:: python

    lvsh3 = A.ls('vsh3')
    lssh3 = A.ls('sh3')
    lmat = A.ls('mat')
    print "Number of antenna in .vsh3 format : ",len(lvsh3)
    print "Number of antenna in .sh3 format : ",len(lssh3)
    print lvsh3[0:5]
    print lssh3[0:5]
    print lmat[0:5]

.. parsed-literal::

    Number of antenna in .vsh3 format :  66
    Number of antenna in .sh3 format :  42
    ['S1R1.vsh3', 'S1R10.vsh3', 'S1R11.vsh3', 'S1R12.vsh3', 'S1R13.vsh3']
    ['S17R1.sh3', 'S17R2m.sh3', 'S1R1.sh3', 'S1R10.sh3', 'S1R11.sh3']
    []


At that point the radiation pattern of the antenna has not yet been
evaluated. The method to calculate the pattern is ``Fsynth()`` with the
``pattern`` option set to true. If the ``pattern`` option is set to
False, the antenna is evaluated for only the specified direction. This
mode is used in the ray tracing technique while the former is used to
visualize the whole antenna pattern.

The vector spherical coefficient are strored in ``A.C``

.. code:: python

    A.C



.. parsed-literal::

    Br
    -------------
    Nf   : 104
    fmin (GHz) : 0.8
    fmax (GHz) : 5.95
    Ncoeff s3 : 72
    
    Bi
    -------------
    Nf   : 104
    fmin (GHz) : 0.8
    fmax (GHz) : 5.95
    Ncoeff s3 : 72
    
    Cr
    -------------
    Nf   : 104
    fmin (GHz) : 0.8
    fmax (GHz) : 5.95
    Ncoeff s3 : 72
    
    Ci
    -------------
    Nf   : 104
    fmin (GHz) : 0.8
    fmax (GHz) : 5.95
    Ncoeff s3 : 72



.. code:: python

    A.Fsynth(pattern=True)
The ``polar()`` method allow to superpose different pattern for a list
of frequencies ``fGHz`` + If ``phd`` (phi in degree) is specified the
diagram is given as a function of :math:`\theta` + If ``thd`` (theta in
degree) is specified the diagram is given as a function of :math:`\phi`

.. code:: python

    f = plt.figure(figsize=(15,15))
    a1 = f.add_subplot(121,polar=True)
    f1,a1 = A.polar(fGHz=[3,4,5],phd=0,GmaxdB=0,fig=f,ax=a1)
    a2 = f.add_subplot(122,polar=True)
    f2,a2 = A.polar(fGHz=[3,4,5],thd=90,GmaxdB=5,fig=f,ax=a2)


.. image:: Antenna_files/Antenna_14_0.png


The vector spherical coefficients can be dispalayed as follows

.. code:: python

    fig = plt.figure(figsize=(8,8))
    A.C.show(typ='s3')
    plt.tight_layout()


.. image:: Antenna_files/Antenna_16_0.png


Defining Antenna gain from analytic formulas
--------------------------------------------

An antenna can also be defined from closed-form expressions. Available
antennas are the following + Omni + Gauss + WirePlate

.. code:: python

    A = Antenna('Omni')
.. code:: python

    A.Fpatt(pattern=True)
.. code:: python

    A.polar()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f2e4ca2ca90>,
     <matplotlib.projections.polar.PolarAxes at 0x7f2e497891d0>)




.. image:: Antenna_files/Antenna_21_1.png


.. code:: python

    A = Antenna('Gauss')
    A.Fsynth()
    A.polar()



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f2e4920f8d0>,
     <matplotlib.projections.polar.PolarAxes at 0x7f2e49209150>)




.. image:: Antenna_files/Antenna_22_1.png

