
PyLayers Overview
=================

Geographical Information
------------------------

-  `Description of Building Layouts <1-GIS/Layout.ipynb>`__
-  `Editing a Layout <1-GIS/LayoutEditor.ipynb>`__
-  `Effects of changing constitutive properties of a
   Layout <1-GIS/Multisubsegments.ipynb>`__
-  `Handling of Geographical Data: The Earth Zone
   class <1-GIS/Ezone.ipynb>`__

Antenna and Propagation
-----------------------

Antennas
~~~~~~~~

-  `Description of Antennas <2-AP/Antenna.ipynb>`__
-  `Scalar Spherical Harmonics Representation <2-AP/AntennaSSH.ipynb>`__
-  `Vector Spherical Harmonics Representation <2-AP/AntennaVSH.ipynb>`__

Slabs and Materials
~~~~~~~~~~~~~~~~~~~

-  `SlabsMaterials <2-AP/SlabsMaterials.ipynb>`__

MultiWall model
~~~~~~~~~~~~~~~

-  `Coverage <2-AP/Coverage.ipynb>`__
-  `MultiWall <2-AP/MultiWall.ipynb>`__

Signatures and Rays
~~~~~~~~~~~~~~~~~~~

-  `The Ray Signatures <2-AP/Signatures.ipynb>`__

-  `The transmission Channel <2-AP/Channel.ipynb>`__

Handling time domain signals
----------------------------

-  `The Bsignal Class <3-PHY/Bsignal.ipynb>`__

Human Mobility
--------------

-  `Simulnet Configuation <5-SIM/SimulNetConfig.ipynb>`__
-  `Large scale mobility ``simulnet`` <4-MOB/Mobility.ipynb>`__

-  `Body scale mobility <4-MOB/Body.ipynb>`__

Ray Tracing Simulation Examples
-------------------------------

-  `Link Simulations <5-SIM/LinkSimulation.ipynb>`__
-  `DLR WHERE2 Scenario <5-SIM/DLR-WHERE2.ipynb>`__
-  `PTIN WHERE2 Scenario <5-SIM/PTIN.ipynb>`__
-  `WHERE1 - M1 Scenario <5-SIM/WHERE1-M1.ipynb>`__
-  `WHERE1 - M1 Scenario (bis) <5-SIM/Where1M1.ipynb>`__
-  `Aggregated CIR <5-SIM/AggregatedCIR.ipynb>`__

PyLayers Classes and Tools
--------------------------

-  `The Cone class <8-MISC/Cone.ipynb>`__
-  `Geometrical Utility Functions <8-MISC/Geomutil.ipynb>`__

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



