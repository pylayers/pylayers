This section explains some features of the class **Layout**. Layout are
handled by the module **pylayers.gis.layout**

In[72]:

.. code:: python

    from pylayers.gis.layout import *
    import  matplotlib.pyplot as plt

Layout ploting
~~~~~~~~~~~~~~


The layout ploting expoits vizaulization function of networkx. The
method showG allows to display the different graphs which constitute the
layout. option graph :math:`$\in$` {'s','t','i','r','v'}

In[73]:

.. code:: python

    L = Layout('TA-Office.str')
    L.dumpr()
    f = plt.figure(figsize=(25,25))
    L.showG(fig=f,graph='s')

Out[73]:

.. parsed-literal::

    (<matplotlib.figure.Figure at 0xac4a02c>,
     <matplotlib.axes.AxesSubplot at 0xa848aec>)

.. image:: 1__Slab_and_Materials_files/1__Slab_and_Materials_fig_00.png

Slab name
---------

Each node of :math:`$\mathcal{G}_s$` with a negative index is a point.

Each node of :math:`$\mathcal{G}_s$` with a positive index corresponds
to a segment (Wall,door,window,...).

Each edge has a **name** attribute, which contains the name entry of its
constitutive Slab.

In[74]:

.. code:: python

    L.Gs.node[20]['name']

Out[74]:

.. parsed-literal::

    'PARTITION'

In the example above the node +20 of :math:`$\mathcal{G}_s$` is a
segment of nature 'PARTITION'.

'PARTITION' is a key entry of the Slab dictionnary.

Slab dictionnary
~~~~~~~~~~~~~~~~

The layout L has a **sl** (slab dictionnary) attribute. Inside are
placed all the available slabs for building the layout.

In[75]:

.. code:: python

    L.sl.keys()

Out[75]:

.. parsed-literal::

    ['WINDOW_GLASS',
     'PLASTERBOARD_7CM',
     'WALL',
     'AIR',
     'WINDOW',
     'METALIC',
     'PLASTERBOARD_14CM',
     'DOOR',
     'FLOOR',
     'METAL',
     'PARTITION',
     'CONCRETE_20CM3D',
     'PLASTERBOARD_10CM',
     'CEIL',
     'CONCRETE_6CM3D',
     'CONCRETE_15CM3D',
     '3D_WINDOW_GLASS',
     'WALLS',
     'WOOD',
     'CONCRETE_7CM3D',
     'PILLAR',
     'ABSORBENT']

Slab Information
~~~~~~~~~~~~~~~~

Each slab contains various information about their constitutive
materials and on electromagnetic propreties.

Below an example for a simple slab, constituted with a single material.
The slab 'WOOD' is made of material 'WOOD'

In[76]:

.. code:: python

    L.sl['WOOD']['lmatname']

Out[76]:

.. parsed-literal::

    ['WOOD']

In[77]:

.. code:: python

    L.sl['WOOD']['thickness']


Out[77]:

.. parsed-literal::

    (4.0, 0, 0, 0, 0, 0, 0, 0)

Until now the thickness attribute is a tuple with 8 values. This is for
backward compatibility with pulsray. In future

pure python version their will remains no reason to maintain this
limitation.

Notice that in thickness the unit is centimeter, lthick is in SI unit
i.e meters.

Ignore the existence of this attribute, use 'lthick' attribute instead.

In[78]:

.. code:: python

    L.sl['WOOD']['lthick']

Out[78]:

.. parsed-literal::

    [0.04]

In[79]:

.. code:: python

    L.sl['WOOD']['color']

Out[79]:

.. parsed-literal::

    'maroon'

In[80]:

.. code:: python

    L.sl['WOOD']['linewidth']

Out[80]:

.. parsed-literal::

    2

More complex Slab, using different stacks of materials can be considered

In[81]:

.. code:: python

    L.sl['3D_WINDOW_GLASS']['lmatname']

Out[81]:

.. parsed-literal::

    ['GLASS', 'AIR', 'GLASS']

In[82]:

.. code:: python

    L.sl['3D_WINDOW_GLASS']['lthick']

Out[82]:

.. parsed-literal::

    [0.005, 0.005, 0.005]

For each material, the electromagnetic propreties can be listed.

In[83]:

.. code:: python

    L.sl['3D_WINDOW_GLASS']['lmat']

Out[83]:

.. parsed-literal::

    [{'epr': (3.79999995232+0j),
      'index': 4,
      'mur': (1+0j),
      'name': 'GLASS',
      'roughness': 0.0,
      'sigma': 0.0},
     {'epr': (1+0j),
      'index': 1,
      'mur': (1+0j),
      'name': 'AIR',
      'roughness': 0.0,
      'sigma': 0.0},
     {'epr': (3.79999995232+0j),
      'index': 4,
      'mur': (1+0j),
      'name': 'GLASS',
      'roughness': 0.0,
      'sigma': 0.0}]

Slab evaluation
~~~~~~~~~~~~~~~


Each Slab can be evaluated to obtain the Transmission and Reflexion
coefficients for

-  a given frequency range
-  a given incidence angle range (:math:`$0\le\theta<\frac{\pi}{2}$`)


In[84]:

.. code:: python

    fGHz   = np.arange(1,5,0.05)
    theta = np.arange(0,pi/2,0.01)
    
    L.sl['WOOD'].ev(fGHz,theta,compensate=True)
    sR = np.shape(L.sl['WOOD'].R) 
    print '\nHere, we have evluated the slab for',sR[0],'frequency(ies)', 'and',sR[1], 'angle(s)\n'

.. parsed-literal::

    
    Here, we have evluated the slab for 80 frequency(ies) and 158 angle(s)
    


Transmission and Reflexion coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the evaluation reflexion and tramsmission coeffcient are computed
for the given frequency range and theta range

In[85]:

.. code:: python

    ifreq=1
    ithet=10
    
    print '\nReflection coefficient @',fGHz[ifreq],'GHz and theta=',theta[ithet],':\n\n R=',L.sl['WOOD'].R[0,0]
    print '\nTransmission coefficient @',fGHz[ifreq],'GHz and theta=',theta[ithet],':\n\n T=',L.sl['WOOD'].T[0,0],'\n'


.. parsed-literal::

    
    Reflection coefficient @ 1.05 GHz and theta= 0.1 :
    
     R= [[-0.46812931-0.06286165j  0.00000000+0.j        ]
     [ 0.00000000+0.j          0.46812931+0.06286165j]]
    
    Transmission coefficient @ 1.05 GHz and theta= 0.1 :
    
     T= [[ 0.72738249-0.48850885j  0.00000000+0.j        ]
     [ 0.00000000+0.j          0.72738249-0.48850885j]] 
    


Ploting Coefficients
~~~~~~~~~~~~~~~~~~~~


 with respect to  frequency
```````````````````````````

In[86]:

.. code:: python

    L.sl['WOOD'].plotwrtf(typ='mod')
    plt.figure()
    L.sl['WOOD'].plotwrtf(typ='phase')

.. image:: 1__Slab_and_Materials_files/1__Slab_and_Materials_fig_01.png

.. image:: 1__Slab_and_Materials_files/1__Slab_and_Materials_fig_02.png

with respect to angle
`````````````````````

In[87]:

.. code:: python

    fGHz= np.array([2.4])
    L.sl['WOOD'].ev(fGHz,theta)
    L.sl['WOOD'].plotwrta()

.. image:: 1__Slab_and_Materials_files/1__Slab_and_Materials_fig_03.png

wrt to angle and frequency 
```````````````````````````

In[91]:

.. code:: python

    figsize(8,8)
    fGHz= np.arange(0.7,5.2,0.1)
    L.sl['WOOD'].ev(fGHz,theta)
    L.sl['WOOD'].pcolor()

.. image:: 1__Slab_and_Materials_files/1__Slab_and_Materials_fig_04.png
