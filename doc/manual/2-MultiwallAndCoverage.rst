Multi-wall model
================


In[10]:

.. code:: python

    import time
    from pylayers.util.project import *
    import pylayers.util.pyutil as pyu
    from pylayers.util.utilnet import str2bool
    from pylayers.gis.layout import Layout
    from pylayers.antprop.multiwall import *
    from pylayers.antprop.coverage import *
    from pylayers.network.model import *


Loading the layout 
-------------------

In[ ]:

.. code:: python

    L=Layout('TA-Office.str')
    try:
        L.dumpr() # load graphs
    except:
        L.build()
        L.dumpw()

Defining a radio link
---------------------

We define :

-  A a transmitter

-  B a receviver



In[11]:

.. code:: python

    A=np.array((4,1)) # transmitter
    B=np.array((30,12)) # receiver

Ploting the scene
-----------------

In[12]:

.. code:: python

    # figure instanciation
    f = plt.figure(figsize=(25,25))
    ax = f.add_subplot(111)
    r = np.array((A,B))
    f,ax = L.showG(fig=f,ax=ax,graph='s')
    
    ax.plot(A[0],A[1],'ob')
    ax.plot(B[0],B[1],'or')
    ax.plot(r[:,0],r[:,1])

Out[12]:

.. parsed-literal::

    [<matplotlib.lines.Line2D at 0xaa60aec>]

.. image:: 2__Multiwall_and_coverage_files/2__Multiwall_and_coverage_fig_00.png

Find the intersection between the "direct" path and the walls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In[13]:

.. code:: python

    L.angleonlink(A,B)

Out[13]:

.. parsed-literal::

    (array([34, 36, 41, 67, 72, 89]),
     array([ 0.40024066,  0.40024066,  0.40024066,  1.17055567,  1.17055567,
            0.40024066]))

Computing the Multi-wall model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The multi-wall model computation returns losses and LOS excess delay for
orthogonal and parallel polarization

In[35]:

.. code:: python

    f   = 2.4
    r= np.array((B,B))
    
    Lwo,Lwp,Edo,Edp=Loss0_v2(L,r,f,A)
    
    print 'Losses polar orthogonal \t %g dBm' %(Lwo[0])
    print 'Losses polar parallel \t %g  dBm' % (Lwp[0])
    print 'Excess delay polar orthogonal \t %g ns' % (Edo[0])
    print 'Excess delay polar parallel  \t %g ns' %Edp[0]

.. parsed-literal::

    Losses polar orthogonal 	 27.7333 dBm
    Losses polar parallel 	 16.0573  dBm
    Excess delay polar orthogonal 	 2.23113 ns
    Excess delay polar parallel  	 2.12364 ns


Coverage class
==============


By extension, the multi-wall model can also been used to perform a
coverage of a Layout

In[15]:

.. code:: python

    C = Coverage()
    C.L  = L # set layout
    C.tx = A # set the tramsitter

In[22]:

.. code:: python

    np.shape(C.grid)

Out[22]:

.. parsed-literal::

    (800, 2)

Compute the coverage
~~~~~~~~~~~~~~~~~~~~


In[16]:

.. code:: python

    t1=time.time()
    C.cover()
    t2=time.time()
    print 'Coverage performed in ', t2-t1, 's'

.. parsed-literal::

    Coverage performed in  5.5845220089 s


Coverage Map
~~~~~~~~~~~~


For Orthogonal polarization

In[17]:

.. code:: python

    C.showPower(polarization='o')
    C.showEd(polarization='o')


.. image:: 2__Multiwall_and_coverage_files/2__Multiwall_and_coverage_fig_01.png

.. image:: 2__Multiwall_and_coverage_files/2__Multiwall_and_coverage_fig_02.png

For parallel polarization

In[18]:

.. code:: python

    C.showPower(polarization='p')
    C.showEd(polarization='p')

.. image:: 2__Multiwall_and_coverage_files/2__Multiwall_and_coverage_fig_03.png

.. image:: 2__Multiwall_and_coverage_files/2__Multiwall_and_coverage_fig_04.png
