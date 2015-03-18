
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



.. code:: python

    from pylayers.signal.bsignal import *
    import numpy as np
.. code:: python

    taumin  = 0
    taumax  = 150
    taustep = 0.1
    x = np.arange(taumin,taumax,taustep)
    y = np.zeros(len(x))
.. code:: python

    CIRa=TUsignal(x,y)
    CIRa



.. parsed-literal::

    TUsignal :  (1500,)  (1500,) 
    time (ns) : 1500



.. code:: python

    CIRa.y



.. parsed-literal::

    array([ 0.,  0.,  0., ...,  0.,  0.,  0.])



.. code:: python

    len(np.shape(CIRa.y))



.. parsed-literal::

    1



.. code:: python

    CIRa.plot(typ='v')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7fc12600bbd0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7fc1260252d0>]], dtype=object))




.. image:: AggregatedCIR_files/AggregatedCIR_6_1.png


.. code:: python

    def aggcir(CIRa,alphak,tauk):
        """ aggregation of CIR from (alphak,tauk)
        
        Parameters
        ----------
        
        alphak : float
            CIR path amplitude
        tauk : float
            CIR delay values
            
        """
        shy = np.shape(CIRa.y)
        x = CIRa.x
        eps = (x[1]-x[0])/2
        u = map(lambda t: np.where( (x>t-eps) & (x<=t+eps))[0][0],tauk)
        ynew  = np.zeros(len(x))
        ynew[u] = alphak
        if len(shy)>1:
           CIRa.y = np.vstack((CIRa.y,ynew))
        else:
           CIRa.y = ynew[np.newaxis,:]
            
        return(CIRa)
.. code:: python

    N = 7
    Ntrial = 100
    for i in range(Ntrial):
        alphak = 10*np.random.rand(N)
        tauk   = taumax *np.random.rand(N)
        CIRa=aggcir(CIRa,alphak,tauk)
.. code:: python

    np.shape(CIRa.y)



.. parsed-literal::

    (100, 1500)



.. code:: python

    plot(CIRa.y[0,:])



.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x7fc125ec6690>]




.. image:: AggregatedCIR_files/AggregatedCIR_10_1.png


.. code:: python

    alphak = 10*np.random.rand(N)
    tauk   = taumax *np.random.rand(N)
    CIRa=aggcir(CIRa,alphak,tauk)
.. code:: python

    np.shape(CIRa.y)



.. parsed-literal::

    (101, 1500)



.. code:: python

    plot(CIRa.y[0,:])



.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x7fc125e12550>]




.. image:: AggregatedCIR_files/AggregatedCIR_13_1.png


.. code:: python

    imshow(CIRa.y)
    axis('tight')



.. parsed-literal::

    (-0.5, 1499.5, 100.5, -0.5)




.. image:: AggregatedCIR_files/AggregatedCIR_14_1.png

