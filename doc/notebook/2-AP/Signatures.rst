
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



Ray Signatures
==============

.. code:: python

    import time
    from pylayers.gis.layout import *
    from pylayers.antprop.signature import *
    from pylayers.antprop.rays import *
    %matplotlib inline

.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


.. code:: python

    L = Layout('defstr3.ini')
    L.build()
    L.dumpw()

::


    ---------------------------------------------------------------------------
    NoOptionError                             Traceback (most recent call last)

    <ipython-input-3-1b28416d0633> in <module>()
    ----> 1 L = Layout('defstr3.ini')
          2 L.build()
          3 L.dumpw()


    /home/uguen/Documents/rch/devel/pylayers/pylayers/gis/layout.pyc in __init__(self, _filename, _filematini, _fileslabini, _filefur, force)
        331         mat.load(_filematini)
        332 
    --> 333         self.sl = sb.SlabDB()
        334         self.sl.mat = mat
        335         self.sl.load(_fileslabini)


    /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/slab.pyc in __init__(self, filemat, fileslab)
       1865             self.mat.load(filemat)
       1866         if (fileslab != ''):
    -> 1867             self.load(fileslab)
       1868             self.dass()
       1869 


    /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/slab.pyc in load(self, _fileini)
       2114         for slabname in self.di.values():
       2115             S=Slab(name=slabname,mat=self.mat)
    -> 2116             S['lmatname']=eval(config.get(slabname,'lmatname'))
       2117             S['nbmat']=len(S['lmatname'])
       2118             S['color']=config.get(slabname,'color')


    /home/uguen/anaconda/lib/python2.7/ConfigParser.pyc in get(self, section, option, raw, vars)
        616             value = d[option]
        617         except KeyError:
    --> 618             raise NoOptionError(option, section)
        619 
        620         if raw or value is None:


    NoOptionError: No option 'lmatname' in section: 'ABSORBENT'


Showing the graph of rooms with 2 rooms separated by a DOOR segment

.. code:: python

    L.showG('sv',figsize=(8,8))
    a=plt.axis('off')


.. image:: Signatures_files/Signatures_5_0.png


The graph of interactions is shown below.

.. code:: python

    L.showG('si',figsize=(20,20))
    a=plt.axis('off')


.. image:: Signatures_files/Signatures_7_0.png


All the interactions of a given cycle are stored as meta information in
nodes of ``Gt``

.. code:: python

    L.Gi.node



.. parsed-literal::

    {(1, 1): {},
     (1, 1, 2): {},
     (1, 2): {},
     (1, 2, 1): {},
     (2, 1): {},
     (2, 1, 2): {},
     (2, 2): {},
     (2, 2, 1): {},
     (3, 1): {},
     (3, 1, 2): {},
     (3, 2): {},
     (3, 2, 1): {},
     (4, 0): {},
     (4, 0, 2): {},
     (4, 2): {},
     (4, 2, 0): {},
     (5, 0): {},
     (5, 0, 2): {},
     (5, 2): {},
     (5, 2, 0): {},
     (6, 0): {},
     (6, 0, 1): {},
     (6, 1): {},
     (6, 1, 0): {},
     (7, 0): {},
     (7, 0, 1): {},
     (7, 1): {},
     (7, 1, 0): {},
     (8, 0): {},
     (8, 0, 1): {},
     (8, 1): {},
     (8, 1, 0): {},
     (9, 0): {},
     (9, 0, 2): {},
     (9, 2): {},
     (9, 2, 0): {}}



.. code:: python

    L.Gt.node[0]['inter']



.. parsed-literal::

    [(6, 0),
     (6, 0, 1),
     (6, 1, 0),
     (7, 0),
     (7, 0, 1),
     (7, 1, 0),
     (8, 0),
     (8, 0, 1),
     (8, 1, 0),
     (9, 0),
     (9, 0, 2),
     (9, 2, 0),
     (4, 0),
     (4, 0, 2),
     (4, 2, 0),
     (5, 0),
     (5, 0, 2),
     (5, 2, 0),
     (-4,),
     (-3,),
     (-1,),
     (-6,)]



The signature is calculated with as parameters the Layout object and two
cycle numbers. In example below it is 0 and 1.

.. code:: python

    Si = Signatures(L,0,1)
The cold start determination of the signature is done with a ``run``
function. The code is not in its final shape here and there is room for
significant acceleration in incorporating propagation based heuristics.
The mitigation of graph exploration depth is done in setting a
``cutoff`` value which limits the exploration in the interaction graph.

.. code:: python

    Si.run5(cutoff=5,diffraction=False,algo='old')
The representation method of a signature gives informations about the
different signatures. Signatures are grouped by number of interactions.

.. code:: python

    L.Gt.pos



.. parsed-literal::

    {0: (758.49, 1111.9),
     1: (766.00300113353387, 1113.947479109665),
     2: (761.00289669547806, 1113.915769812613)}



.. code:: python

    ptx = np.array(L.Gt.pos[0])+np.random.rand(2)
    prx = np.array(L.Gt.pos[1])+np.random.rand(2)
    print ptx
    print prx

.. parsed-literal::

    [  759.31624319  1112.04712068]
    [  766.78504482  1114.60018036]

