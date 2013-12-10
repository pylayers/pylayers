
Ray Signatures
==============


.. code-block:: python

    import time
    from pylayers.gis.layout import *
    from pylayers.antprop.signature import *
    from pylayers.antprop.rays import *


.. parsed-literal::

    <matplotlib.figure.Figure at 0x39b0050>


.. code-block:: python

    L = Layout()
    try:
        L.dumpr()
    except:
        L.build()
        L.dumpw()
Showing the graph of rooms with 2 rooms separated by a DOOR segment

.. code-block:: python

    L.showG('v')
    a=plt.axis('off')


.. image:: Signatures_files/Signatures_4_0.png


The graph of interactions is shown below.

.. code-block:: python

    L.showG('i',figsize=(20,20))
    a=plt.axis('off')


.. image:: Signatures_files/Signatures_6_0.png


All the interactions of a given cycle are stored as meta information in
nodes of ``Gt``

.. code-block:: python

    L.Gt.node[0]['inter']



.. parsed-literal::

    ['(3, 0)',
     '(4, 0)',
     '(7, 0)',
     '(7, 0, 1)',
     '(7, 1, 0)',
     '(9, 0)',
     '(9, 0, 1)',
     '(9, 1, 0)',
     '(8, 0)',
     '(8, 0, 1)',
     '(8, 1, 0)',
     '(2, 0)']



The signature is calculated with as parameters the Layout object and two
cycle numbers. In example below it is 0 and 1.

.. code-block:: python

    Si = Signatures(L,0,1)
The cold start determination of the signature is done with ``run1``
function. The code is not finished here and there is room for
significant acceleration in incorporating propagation based heuristics.
The mitigation of combinatorics explosion is done in setting a cutoff
vlaue which limits the exploration in the interaction graph.

.. code-block:: python

    Si.run1(cutoff=3)
An exhaustive search of signatures when no prior information is given is
a combinatory hard problem. The number of signatures grows rapidly with
the cutoff parameter.

The representaion method of a signature gives informations about the
different signatures. Signatures are grouped by number of interactions.

.. code-block:: python

    Si



.. parsed-literal::

    Signatures
    ----------
    from cycle : 0 to cycle 1
    1 : 3
       [7 9 8]
       [2 2 2]
    2 : 18
       [7 7 7 9 9 9 8 8 8 3 3 3 4 4 4 2 2 2]
       [2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1]
       [6 5 1 6 5 1 6 5 1 7 9 8 7 9 8 7 9 8]
       [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2]
    3 : 69
       [7 7 7 7 7 7 7 9 9 9 9 9 9 9 9 9 8 8 8 8 8 8 8 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     3 4 4 4 4 4 4 4 4 4 4 4 4 7 7 9 9 9 8 8 2 2 2 2 2 2 2 2 2 2 2 2]
       [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
       [5 1 1 6 6 5 6 5 1 6 1 6 6 6 5 6 5 1 6 1 6 6 5 4 2 4 2 4 2 9 8 7 9 8 7 9 8
     7 2 3 2 3 2 9 7 8 7 9 7 8 3 3 3 3 3 3 3 3 4 3 4 4 9 8 7 9 8 7 8]
       [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
     2 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2]
       [6 6 5 9 8 1 1 6 6 5 5 7 9 8 1 1 6 6 5 5 7 9 1 7 7 9 9 8 8 6 6 6 5 5 5 1 1
     1 7 9 9 8 8 6 6 6 5 1 1 1 9 8 7 9 8 7 9 7 7 9 9 8 6 6 6 5 5 5 1]
       [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1
     1 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1]
    4 : 258
       [7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
     9 9 9 9 9 9 9 9 9 9 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 3 3 3 3 3 3
     3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
     4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 9 9 9 9 9 9 9 9 9 9 9 9 9
     9 9 9 9 9 9 9 9 9 9 9 9 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 2 2 2 2 2 2
     2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
       [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
       [5 1 6 6 1 5 6 1 1 5 5 6 1 1 5 5 6 1 5 6 6 5 1 6 6 6 1 5 6 6 1 1 5 5 6 6 1
     1 5 5 6 6 1 5 6 6 6 5 1 6 6 1 5 6 1 1 5 5 6 1 1 5 5 6 1 5 6 6 4 2 4 2 4 2
     9 9 8 8 4 4 4 7 7 2 2 2 9 9 8 8 4 7 2 2 2 9 8 9 8 7 9 7 9 9 8 4 4 4 7 7 2
     3 2 9 3 7 2 9 3 7 9 3 3 9 7 7 2 2 2 8 3 3 9 9 7 2 2 2 8 8 9 8 9 7 8 9 7 3
     3 9 7 7 2 3 2 2 4 3 2 2 4 4 3 2 4 4 3 3 3 3 3 3 3 3 2 2 4 3 3 2 2 4 4 3 3
     2 4 4 3 3 3 3 3 3 3 3 3 3 2 2 4 3 2 2 4 4 3 2 4 4 3 3 3 3 3 3 3 8 9 3 4 8
     9 3 4 9 3 3 9 8 8 4 4 4 7 3 3 9 8 8 4 9 8 9 8 7 9 7 3 3 9 9 8 4 4 4 7 7]
       [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     2 2 2 2 1 1 1 2 2 1 1 1 2 2 2 2 1 2 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 1
     1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 2 1 1 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 1
     1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 2 2 2 1 1 1 2 1 1 2 2 2 1 2 2 2 2 2 2 2 1 1 2 2 2 1 1 1 2 2]
       [1 5 9 1 6 1 1 5 6 1 6 1 5 6 1 6 1 5 6 9 8 1 5 7 9 1 6 1 1 5 5 6 1 6 1 5 5
     6 1 6 1 5 5 6 9 8 5 1 5 7 9 6 1 5 5 6 1 6 5 5 6 1 6 5 5 6 9 5 2 4 2 4 2 4
     5 1 5 1 9 7 8 5 1 9 8 7 6 1 6 1 7 1 9 8 7 6 6 6 6 6 6 6 5 6 5 9 7 8 5 6 8
     2 3 3 2 3 3 3 2 3 3 9 8 1 5 1 9 8 7 1 9 8 6 1 1 9 8 7 6 1 6 6 6 6 6 6 6 9
     8 6 5 6 8 2 3 4 2 2 3 4 3 2 2 4 3 2 9 8 9 8 9 8 4 2 3 4 2 4 2 3 4 3 2 4 2
     4 3 2 9 8 7 9 8 7 9 8 7 4 3 4 2 4 3 4 3 2 4 4 3 2 9 7 9 7 9 7 4 3 3 4 3 3
     3 4 3 3 9 7 5 5 1 9 7 8 5 9 7 6 6 1 7 6 6 6 6 6 6 6 9 7 5 6 5 9 7 8 5 6]
       [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 2 2 2 1 1 2 2 2 1 1 1 1 2 1 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 2
     1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 2 2 2 1 2 2 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 2
     2 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
     1 1 1 1 2 2 1 1 1 2 2 2 1 2 2 1 1 1 2 1 1 1 1 1 1 1 2 2 1 1 1 2 2 2 1 1]
       [6 6 5 5 5 7 7 7 7 9 9 9 9 9 8 8 8 8 1 1 1 6 6 5 5 5 5 7 7 7 7 7 9 9 9 9 9
     9 8 8 8 8 8 1 1 1 1 6 6 5 5 5 7 7 7 7 9 9 9 9 9 8 8 8 8 1 1 1 7 7 9 9 8 8
     6 6 6 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 7 7 9 9 9 8 8 1 1 1 1 1 1 1 1 1
     7 7 7 9 9 9 9 8 8 8 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 7 7 9 9 9 8 8 1
     1 1 1 1 1 7 7 7 7 9 9 9 9 9 8 8 8 8 6 6 5 5 1 1 7 7 7 7 7 9 9 9 9 9 9 8 8
     8 8 8 6 6 6 5 5 5 1 1 1 7 7 7 7 9 9 9 9 9 8 8 8 8 6 6 5 5 1 1 7 7 7 9 9 9
     9 8 8 8 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 7 7 9 9 9 8 8 1 1 1 1 1 1 1 1 1 1]
       [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2
     2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2
     2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]




.. code-block:: python

    L.Gt.pos



.. parsed-literal::

    {0: (7.5, 0.0), 1: (2.5, -0)}



.. code-block:: python

    ptx = np.array(L.Gt.pos[0])+np.random.rand(2)
    prx = np.array(L.Gt.pos[1])+np.random.rand(2)
    print ptx
    print prx

.. parsed-literal::

    [ 7.88647153  0.322159  ]
    [ 2.88164786  0.68333672]


Evaluate performances of signature and ray evaluation

.. code-block:: python

    tt1 = {}
    tt2 = {}
    tint={}
    tsig={} # number of sig
    tray={} # number of rays
    maxcutoff=7
    for cutoff in range(maxcutoff):
         del Si
            
         Si = Signatures(L,0,1)
         tic1=time.time()
         Si.run1(cutoff=cutoff)
         tic2=time.time()
         r2d = Si.rays(ptx,prx)   
         tic3=time.time()
         Si.num()   
         Nr = len(r2d)   
         tt1[cutoff]=tic2-tic1
         tt2[cutoff]=tic3-tic2   
         tint[cutoff]=Si.nint   
         tsig[cutoff]=Si.nsig
         tray[cutoff]=Nr
         print cutoff,tt2[cutoff]
    #    L.display['ednodes']=False
    #    r2d.show(L)

.. parsed-literal::

    0 0.00191807746887
    1 0.0103960037231
    2 0.0516300201416
    3 0.208976984024
    4 0.682647943497
    5 1.92105197906
    6 4.61728715897


It appears that the increasing number of signatures obatined with the
run1 algorithm do not yield necessarily a significant increase in the
number of rays.

.. code-block:: python

    p1=semilogy(tt1.keys(),tray.values(),'ob')
    p2=semilogy(tt1.keys(),tsig.values(),'or')
    legend((p1[0],p2[0]),('rays','signature'),loc='best')
    xlabel('cutoff')
    ylabel('#')



.. parsed-literal::

    <matplotlib.text.Text at 0x4817790>




.. image:: Signatures_files/Signatures_21_1.png


.. code-block:: python

    #b1=bar(tt2.keys(),tt2.values(),color='red')
    #b2=bar(tt1.keys(),tt1.values(),color='blue')
    #b1=semilogx(tt2.keys(),log10(tt2.values()),color='red')
    #b2=semilogx(tt1.keys(),log10(tt1.values()),color='blue')
    #b3=semilogx(tt1.keys(),log10(0.1*(arange(10))**3),color='green')
    b1=loglog(tt2.keys(),tt2.values(),'k.')
    b2=loglog(tt1.keys(),tt1.values(),'r.')
    b3=loglog(tt1.keys(),0.1*arange(maxcutoff)**3,color='green')
    b4=loglog(tt1.keys(),0.001*arange(maxcutoff)**2,color='green')
    b5=loglog(tt1.keys(),0.1*arange(maxcutoff)**4.4,color='cyan')
    legend((b1[0],b2[0],b3[0],b4[0],b5[0]),('rays','signature',r'$O(N^3)$',r'$O(N^2)$',r'$O(N^{4.4})$'),loc='best')
    xlabel('cutoff')
    ylabel('time (s)')



.. parsed-literal::

    <matplotlib.text.Text at 0x4a081d0>




.. image:: Signatures_files/Signatures_22_1.png


.. code-block:: python

    b1=loglog(tt2.keys(),tint.values(),color='red')
    b2=loglog(tt1.keys(),tsig.values(),color='blue')
    #b3=loglog(tt1.keys(),arange(maxcutoff)**3,color='green')
    #b4=loglog(tt1.keys(),arange(maxcutoff)**2,color='green')
    #legend((b1[0],b2[0],b3[0],b4[0]),('interactions','signatures',r'$O(N^3)$',r'$O(N^2)$'),loc='best')
    legend((b1[0],b2[0]),('interactions','signatures'),loc='best')
    xlabel('cutoff')
    ylabel('#')



.. parsed-literal::

    <matplotlib.text.Text at 0x4b1aad0>




.. image:: Signatures_files/Signatures_23_1.png


.. code-block:: python

    #b1=loglog(tsig.values(),tt1.values(),color='red')
    b2=loglog(tint.values(),tt2.values(),color='red')
    b3=loglog(tint.values(),0.0004*array(tint.values())**0.88,color='blue')
    xlabel('# Interactions')



.. parsed-literal::

    <matplotlib.text.Text at 0x7ff7f4c066d0>




.. image:: Signatures_files/Signatures_24_1.png


.. code-block:: python

    plot(tint.values(),tt2.values(),color='red')
    plot(tint.values(),0.0004*array(tint.values())**0.88,color='blue')



.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x7ff7f4f9fa50>]




.. image:: Signatures_files/Signatures_25_1.png


The computation time in rays grows exponentially with the number of
interactions.

:math:`$$T_{rays}=\alpha N_i^{0.88} [s]$$`

.. code-block:: python

    b1=loglog(tt2.keys(),tint.values(),color='red')
    b2=loglog(tt2.keys(),array(tt2.keys())**5,color='blue')


.. image:: Signatures_files/Signatures_27_0.png


.. code-block:: python

    for k in range(maxcutoff):
        r2d.show(L,i=k+1,colray='red',widthray=0.5)
        title(str(k+1))


.. image:: Signatures_files/Signatures_28_0.png



.. image:: Signatures_files/Signatures_28_1.png



.. image:: Signatures_files/Signatures_28_2.png



.. image:: Signatures_files/Signatures_28_3.png



.. image:: Signatures_files/Signatures_28_4.png



.. image:: Signatures_files/Signatures_28_5.png



.. image:: Signatures_files/Signatures_28_6.png


.. code-block:: python

    fig,ax=r2d.show(L,i=1,figsize=(20,10),colray='red',widthray=3)
    fig,ax=r2d.show(L,i=1,colray='green',widthray=2,fig=fig,ax=ax)
    fig,ax=r2d.show(L,i=2,colray='blue',widthray=1,fig=fig,ax=ax)
    fig,ax=r2d.show(L,i=3,colray='black',widthray=1,fig=fig,ax=ax)
    fig,ax=r2d.show(L,i=4,colray='green',widthray=0.5,fig=fig,ax=ax)
    fig,ax=r2d.show(L,i=5,colray='green',widthray=0.5,fig=fig,ax=ax)
    fig,ax=r2d.show(L,i=6,colray='green',widthray=0.5,fig=fig,ax=ax)
    fig,ax=r2d.show(L,i=7,colray='green',widthray=0.5,fig=fig,ax=ax)
    #fig,ax=r2d.show(L,i=8,colray='green',widthray=0.5,fig=fig,ax=ax)
    #fig,ax=r2d.show(L,i=9,colray='green',widthray=0.5,fig=fig,ax=ax)


.. image:: Signatures_files/Signatures_29_0.png


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


