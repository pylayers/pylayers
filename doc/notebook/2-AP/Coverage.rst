
Indoor Coverage with the PyLayers Multi Wall model
==================================================

    A Multi wall model accounts only for the attenuation along the
    direct path between Tx and Rx.

    This approach do not provide any information about delay spread or
    multi paths but it can nevertheless be useful in different context
    as optimization of an indoor radio network. The MultiWall approach
    provides a fast indication about the propagation losses.

    A ray tracing approach is much more precise, but also is much more
    time consuming and depending on the purpose, it could be relevant to
    proceed with a simpler and faster site-specific approach as the
    Multiwall model.

    **``PyLayers``** provides a multiwall module which heavily relies on
    the core class **``Slab``**. Notice that, the same core **``Slab``**
    module is used for Ray tracing and MultiWall model approaches.

Let see, how it works. First let's import the ``coverage`` module. And
the ``time`` module for performance evaluation.

.. code:: python

    from pylayers.antprop.coverage import *
    import time

.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


Instantiate a coverage object. By defaut the ``TA-Office.ini`` layout is
loaded.

The coverage information is placed in the **coverage.ini** file in the
project directory.

Below is an example such a configuration file.

.. code:: python

    !cat $BASENAME/ini/coverage.ini

.. parsed-literal::

    [grid]
    nx = 40
    ny = 20
    boundary = [20,0,30,20]
    mode = full ; file, zone , full 
    file = 'points.ini'
    
    [layout]
    filename = TA-Office.ini ; 0 40 0 15
    ;filename = W2PTIN.ini
    ;filename = Lstruc.str
    
    [ap]
    0 = {'name':'room1','wstd':'ieee80211b','p':(1,12,1.2),'PtdBm':0,'chan':[11],'on':True} 
    1 = {'name':'room2','wstd':'ieee80211b','p':(10,2,1.2),'PtdBm':0,'chan':[11],'on':True} 
    2 = {'name':'room3','wstd':'ieee80211b','p':(20,1,1.2),'PtdBm':0,'chan':[11],'on':True} 
    3 = {'name':'room4','wstd':'ieee80211b','p':(36.5,1.5,1.2),'PtdBm':0,'chan':[11],'on':True} 
    4 = {'name':'room5','wstd':'ieee80211b','p':(25,12,1.2),'PtdBm':0,'chan':[11],'on':True} 
    
    [rx]
    temperaturek = 300
    noisefactordb = 0 
    
    [show]
    show = True


The ini file contains 5 sections.

-  [grid] section

   This section precises the size of the grid. By default the grid is
   placed over the whole region of the Layout. A selected region can
   also be defined whith the ``boundary`` list

-  [layout] section

   The name of the layout file (filename = )

-  [ap] section

   A dictionnary of access points precising the standard, the used
   channel, the emitted power and the position of the access point.

-  [show] section

.. code:: python

    # Create a Coverage object from coverag.ini file
    C = Coverage('coverage.ini')
``C`` has a dictionnary ``dap`` (dictionnary of access points) which
gathers information about each access points of the scene.

.. code:: python

    C.dap[1]



.. parsed-literal::

    name : room2
    p : (10, 2, 1.2)
    PtdBm : 0
    channels  : [11]   2.462 : [2.451,2.473]
    sensdBm : -94
    nant : 1
    On : True




The coverage object has a ``__repr__`` method which summarizes different
parameters of the current coverage object

.. code:: python

    C



.. parsed-literal::

    Layout file : TA-Office.ini
    
    -----list of Access Points ------
    name : room1
    p : (1, 12, 1.2)
    PtdBm : 0
    channels  : [11]   2.462 : [2.451,2.473]
    sensdBm : -94
    nant : 1
    On : True
    
    name : room2
    p : (10, 2, 1.2)
    PtdBm : 0
    channels  : [11]   2.462 : [2.451,2.473]
    sensdBm : -94
    nant : 1
    On : True
    
    name : room3
    p : (20, 1, 1.2)
    PtdBm : 0
    channels  : [11]   2.462 : [2.451,2.473]
    sensdBm : -94
    nant : 1
    On : True
    
    name : room4
    p : (36.5, 1.5, 1.2)
    PtdBm : 0
    channels  : [11]   2.462 : [2.451,2.473]
    sensdBm : -94
    nant : 1
    On : True
    
    name : room5
    p : (25, 12, 1.2)
    PtdBm : 0
    channels  : [11]   2.462 : [2.451,2.473]
    sensdBm : -94
    nant : 1
    On : True
    
    -----Rx------
    temperature (K) : 300
    noisefactor (dB) : 0
    
    --- Grid ----
    mode : full
    nx : 40
    ny : 20



Then, the coverage calculation is launched by calling the ``cover()``
method

.. code:: python

    tic = time.time()
    C.cover()
    toc = time.time()
    print "Execution time : %2.3f " % (toc-tic) 

.. parsed-literal::

    Execution time : 2.880 


Let display the current Layout with hidding nodes.

.. code:: python

    from matplotlib.pyplot import *
    %matplotlib inline
    fig=figure(figsize=(10,5))
    C.L.display['nodes']=False
    C.L.display['ednodes']=False
    f,a = C.show(fig=fig)


.. image:: Coverage_files/Coverage_15_0.png


The shadowing map coverage results can be displayed by invoquing various
functions.

.. code:: python

    fig=figure(figsize=(10,5))
    f,a=C.show(fig=fig,typ='pr')


.. image:: Coverage_files/Coverage_17_0.png


.. code:: python

    fig=figure(figsize=(10,5))
    f,a=C.show(fig=fig,typ='pr',f=4)


.. image:: Coverage_files/Coverage_18_0.png


.. code:: python

    fig=figure(figsize=(10,5))
    f,a=C.show(fig=fig,typ='pr',f=10)


.. image:: Coverage_files/Coverage_19_0.png


.. code:: python

    fig=figure(figsize=(10,5))
    f,a=C.show(fig=fig,typ='best',f=1)


.. image:: Coverage_files/Coverage_20_0.png


.. code:: python

    fig=figure(figsize=(10,5))
    f,a=C.show(fig=fig,typ='best',f=10)


.. image:: Coverage_files/Coverage_21_0.png


.. code:: python

    fig=figure(figsize=(10,5))
    C.show(fig=fig,f=5,typ='sinr')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f58c9a66a50>,
     <matplotlib.axes.AxesSubplot at 0x7f58cb7a4710>)




.. image:: Coverage_files/Coverage_22_1.png


As you have noticed the calculation has been done for all the center
frequencies of the selected standard. This is done in prevision of
further channel optimizations.

Let's consider an other standard

.. code:: python

    C2 = Coverage('coverage2.ini')
    C2.cover()
.. code:: python

    C2.show(ftyp='pr')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f58c9a0fa50>,
     <matplotlib.axes.AxesSubplot at 0x7f58c9459590>)




.. image:: Coverage_files/Coverage_26_1.png


.. code:: python

    C.snro.shape



.. parsed-literal::

    (13, 800, 5)



.. code:: python

    fig=figure(figsize=(10,5))
    C.show(fig=fig,f=5,typ='capacity',dB=False)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f58c95295d0>,
     <matplotlib.axes.AxesSubplot at 0x7f58c9529d90>)




.. image:: Coverage_files/Coverage_28_1.png


All simulated quantities are stored in linear scale.

.. code:: python

    C2.Lwo[0,0,0]



.. parsed-literal::

    0.078045027166146197



.. code:: python

    C2.freespace[0,0,0]



.. parsed-literal::

    6.7682907399583888e-07


