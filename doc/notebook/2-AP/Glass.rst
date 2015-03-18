
.. code:: python

    %pylab inline
    from pylayers.antprop.slab import *

.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib


Slabs and Materials
===================



A slab is a set ol several layers of materials with specified thickness.
Slabs are used to describe properties of the different constitutive
elements of a building such as wall, windows ,...

In practice when describing a specific building, it is necessary to
specify a set of different slabs with different characteristics.

The structure which gathers this set is ``SlabDB``. If no file argument
is given, this structure is initialized with the default file:
`slabDB.ini <https://github.com/pylayers/pylayers/blob/master/data/ini/slabDB.ini>`__

This section demonstrates some features of the ``pylayers.antprop.slab``
module.

The Class ``SlabDB`` contains a dictionnary of all available Slabs. This
information is read in the file ``slabDB.ini`` of the current project
pointed by environment variable ``$BASENAME``

.. code:: python

    S = SlabDB()
A SlabDB is a dictionnary, the keys are for the current file are shown
below

Defining a new Slab and a new Material
--------------------------------------

.. code:: python

    S.mat.add(name='glass',typ='epsr',cval=6.76,sigma=0)
    S.mat.add(name='Argent',typ='epsr',cval=1,sigma=63e6)
`Modele de Drude <http://fr.wikipedia.org/wiki/Mod%C3%A8le_de_Drude>`__

.. code:: python

    S.add('s1',['glass'],[0.006])
    S.add('s1Ag',['glass','Argent'],[0.006,7e-9])
    S.add('s2Ag',['glass','Argent','AIR','Argent'],[0.006,4e-9,1e-10,4e-9])
    S.add('s3Ag',['glass','Argent','AIR','Argent','AIR','Argent'],[0.006,4e-9,1e-10,4e-9,1e-10,4e-9])
    S.add('DGs1',['glass','AIR','glass'],[0.008,0.012,0.006])
    S.add('DGs12Ag',['glass','AIR','Argent','AIR','Argent','glass'],[0.008,0.012,4e-9,1e-10,4e-9,0.006])
    S.add('TGs1',['glass','AIR','glass','AIR','glass'],[0.004,0.012,0.004,0.012,0.004])
    S.add('TGs12Ag1',['glass','Argent','AIR','glass','AIR','Argent','glass'],[0.004,7e-9,0.012,0.004,0.012,7e-9,0.004])
.. code:: python

    S.mat['Argent']



.. parsed-literal::

    {'epr': 1,
     'fGHz': 1,
     'index': 12,
     'mur': 1,
     'name': 'Argent',
     'roughness': 0,
     'sigma': 63000000.0}



.. code:: python

    S['s1Ag']['lmatname']



.. parsed-literal::

    ['glass', 'Argent']



.. code:: python

    S['s1Ag']['lthick']



.. parsed-literal::

    [0.006, 7e-09]



.. code:: python

    fGHz= np.arange(0.5,18,0.3)
    theta = np.arange(0,np.pi/2,0.01)
    S['s1'].ev(fGHz,theta)
    S['s1Ag'].ev(fGHz,theta)
    S['s2Ag'].ev(fGHz,theta)
    S['s3Ag'].ev(fGHz,theta)
    S['DGs1'].ev(fGHz,theta)
    S['DGs12Ag'].ev(fGHz,theta)
    S['TGs1'].ev(fGHz,theta)
    S['TGs12Ag1'].ev(fGHz,theta)
.. code:: python

    #S['s1'].plotwrt(var='f',coeff='T',polar='o')
    S['s1'].plotwrt(var='f',coeff='T',polar='po',att=False)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f143d3250>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f143d36d0>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f142ed9d0>]], dtype=object))




.. image:: Glass_files/Glass_15_1.png


.. code:: python

    U =  S['s1'].T[:,:,1,1]
.. code:: python

    U.shape



.. parsed-literal::

    (59, 158)



.. code:: python

    f=S['s1'].fGHz
.. code:: python

    W=U.T
.. code:: python

    from pylayers.signal.bsignal import *
.. code:: python

    T=FUsignal(f,W)
.. code:: python

    fig=figure(figsize(10,10))
    T.show(fig=fig,fontsize=12,cmap=cm.hot,vmin=-15)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f05340290>,
     [<matplotlib.axes.AxesSubplot at 0x7f9f053402d0>,
      <matplotlib.axes.AxesSubplot at 0x7f9f0521acd0>])




.. image:: Glass_files/Glass_22_1.png


.. code:: python

    S['s1Ag'].plotwrt(var='f',coeff='R',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f0513ff50>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f051ba1d0>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f05074d50>]], dtype=object))




.. image:: Glass_files/Glass_23_1.png


.. code:: python

    S['s1Ag'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f050c61d0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f051c2950>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f04d96590>]], dtype=object))




.. image:: Glass_files/Glass_24_1.png


.. code:: python

    S['s2Ag'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f04df5dd0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f04ba63d0>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f04b45a90>]], dtype=object))




.. image:: Glass_files/Glass_25_1.png


.. code:: python

    S['s3Ag'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f04a57790>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f049a4b50>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f04a2fe90>]], dtype=object))




.. image:: Glass_files/Glass_26_1.png


.. code:: python

    S['DGs1'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f048a3cd0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f0495d510>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f0486f2d0>]], dtype=object))




.. image:: Glass_files/Glass_27_1.png


.. code:: python

    S['DGs12Ag'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f04707250>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f046e3d90>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f046727d0>]], dtype=object))




.. image:: Glass_files/Glass_28_1.png


.. code:: python

    S['TGs1'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f0452f4d0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f048a3e10>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f044c5bd0>]], dtype=object))




.. image:: Glass_files/Glass_29_1.png


.. code:: python

    S['TGs12Ag1'].plotwrt(var='f',coeff='T',polar='po',att=True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f9f0430dc90>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f9f04313290>],
            [<matplotlib.axes.AxesSubplot object at 0x7f9f042f2150>]], dtype=object))




.. image:: Glass_files/Glass_30_1.png


.. code:: python

    figure(figsize=(10,10))
    S['TGs12Ag1'].pcolor()


.. image:: Glass_files/Glass_31_0.png


.. code:: python

    figure(figsize=(10,10))
    S['TGs1'].pcolor()


.. image:: Glass_files/Glass_32_0.png


As any PyLayers object there is an help function for remembering which
methods are implemented in the class.

.. code:: python

    S['WOOD']['lmatname']



.. parsed-literal::

    ['WOOD']



thickness is expressed in meters

.. code:: python

    S['WOOD']['lthick']



.. parsed-literal::

    [0.04]



.. code:: python

    S['WOOD']['color']



.. parsed-literal::

    'maroon'



.. code:: python

    S['WOOD']['linewidth']



.. parsed-literal::

    2



Multi layers Slab, using different stacks of materials can be easily
defined using the two lists **lmatname** and **lthick**.

    Notice the adopted convention naming lists starting with letter 'l'
    and dictionnaries starting with letter 'd'

.. code:: python

    S['3D_WINDOW_GLASS']['lmatname']



.. parsed-literal::

    ['GLASS', 'AIR', 'GLASS']



.. code:: python

    S['3D_WINDOW_GLASS']['lthick']



.. parsed-literal::

    [0.005, 0.005, 0.005]



For each constitutive material of a slab, their electromagnetic
properties can be obtained as:

.. code:: python

    S['3D_WINDOW_GLASS']['lmat']



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



Evaluation of a Slab
--------------------

Each Slab can be evaluated to obtain the Transmission and Reflexion
coefficients for

-  a given frequency range
-  a given incidence angle range (:math:`0\le\theta<\frac{\pi}{2}`)

.. code:: python

    fGHz = np.arange(3,5,0.01)
    theta = np.arange(0,np.pi/2,0.01)
    
    S['WOOD'].ev(fGHz,theta,compensate=True)
    sR = np.shape(S['WOOD'].R) 
    print '\nHere, slab is evaluted for',sR[0],'frequency(ies)', 'and',sR[1], 'angle(s)\n'

.. parsed-literal::

    
    Here, slab is evaluted for 200 frequency(ies) and 158 angle(s)
    


Transmission and Reflexion coefficients
---------------------------------------

Reflexion and transmission coefficient are computed for the given
frequency range and theta range

.. code:: python

    ifreq=1
    ithet=10
    
    print '\nReflection coefficient @',fGHz[ifreq],'GHz and theta=',theta[ithet],':\n\n R=',S['WOOD'].R[0,0]
    print '\nTransmission coefficient @',fGHz[ifreq],'GHz and theta=',theta[ithet],':\n\n T=',S['WOOD'].T[0,0],'\n'


.. parsed-literal::

    
    Reflection coefficient @ 3.01 GHz and theta= 0.1 :
    
     R= [[-0.39396205-0.17289585j  0.00000000+0.j        ]
     [ 0.00000000+0.j          0.39396205+0.17289585j]]
    
    Transmission coefficient @ 3.01 GHz and theta= 0.1 :
    
     T= [[-0.17594898-0.86927604j -0.00000000+0.j        ]
     [-0.00000000+0.j         -0.17594898-0.86927604j]] 
    


Ploting Reflection and Transmission Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The method ``plotwrt`` can plot the different calculated coefficients
with respect to angle or frequency.

.. code:: python

    S['WOOD']['lthick']=[0.02]
    S['WOOD'].ev()
    S['WOOD'].ev()
    f,a=S['WOOD'].plotwrt()


.. image:: Glass_files/Glass_52_0.png


.. code:: python

    fGHz = np.arange(1,10,0.01)
    theta = np.arange(0,np.pi/2,0.01)
    
    S['3D_WINDOW_GLASS']['lthick']=[0.006,0.01,0.006]
    #S['3D_WINDOW_GLASS']['lmatname']=['GLASS','AIR','GLASS']
    S['3D_WINDOW_GLASS'].ev(fGHz,theta)
.. code:: python

    fig,ax = S['3D_WINDOW_GLASS'].plotwrt(var='f',coeff='T',polar='o')


.. image:: Glass_files/Glass_54_0.png


.. code:: python

    fig,ax = S['WOOD'].plotwrt(var='a',coeff='R',polar='p')


.. image:: Glass_files/Glass_55_0.png


plot with respect to angle

.. code:: python

    fig = plt.figure(figsize=(20,20))
    fGHz= np.array([2.4])
    S['WOOD'].ev(fGHz,theta)
    fig,ax = S['WOOD'].plotwrt(var='a',coeff='R',fig=fig)
    plt.tight_layout()



.. parsed-literal::

    <matplotlib.figure.Figure at 0x7f9f029b8790>



.. image:: Glass_files/Glass_57_1.png


wrt to angle and frequency

.. code:: python

    plt.figure(figsize=(10,10))
    fGHz= np.arange(0.7,5.2,0.1)
    S['WOOD'].ev(fGHz,theta)
    S['WOOD'].pcolor()


.. image:: Glass_files/Glass_59_0.png


.. code:: python

    theta = np.arange(0,np.pi/2,0.01)
    fGHz = np.arange(0.1,10,0.2)
    sl = SlabDB('matDB.ini','slabDB.ini')
    mat   = sl.mat
    lmat  = [mat['AIR'],mat['WOOD']]
    II    = MatInterface(lmat,0,fGHz,theta)
    II.RT()
    fig,ax = II.plotwrt(var='a',kv=10,typ=['m'])
    plt.tight_layout()
    air = mat['AIR']
    brick  = mat['BRICK']
    II  = MatInterface([air,brick],0,fGHz,theta)
    II.RT()
    fig,ax = II.plotwrt(var='f',color='k',typ=['m'])
    plt.tight_layout()


.. image:: Glass_files/Glass_60_0.png



.. image:: Glass_files/Glass_60_1.png


.. code:: python

    ## Adding new materials
.. code:: python

    sl.mat.add(name='TESS-p50',cval=3+0j,sigma=0.06,typ='epsr')
    
    sl.add(name='TESS-p50-5cm',lmatname=['TESS-p50'],lthick=[0.05])
    sl.add(name='TESS-p50-10cm',lmatname=['TESS-p50'],lthick=[0.10])
    sl.add(name='TESS-p50-15cm',lmatname=['TESS-p50'],lthick=[0.15])
    fGHz=4
    theta = np.arange(0,np.pi/2,0.01)
    #figure(figsize=(8,8))
    # These Tessereau page 50 
    
    sl['TESS-p50-5cm'].ev(fGHz,theta,compensate=True)
    sl['TESS-p50-10cm'].ev(fGHz,theta,compensate=True)
    sl['TESS-p50-15cm'].ev(fGHz,theta,compensate=True)
    
    # by default var='a' and kv = 0 
    
    fig,ax = sl['TESS-p50-5cm'].plotwrt(color='k',labels=['5cm'])
    fig,ax = sl['TESS-p50-10cm'].plotwrt(color='k',labels=['10cm'],linestyle='dashed',fig=fig,ax=ax)
    fig,ax = sl['TESS-p50-15cm'].plotwrt(color='k',labels=['15cm'],linestyle='dashdot',fig=fig,ax=ax)
    plt.axis([0,90,-40,0])
    plt.tight_layout()


.. image:: Glass_files/Glass_62_0.png


Evaluation without phase compensation
-------------------------------------

.. code:: python

    fGHz = np.arange(2,16,0.1)
    theta = 0 
    
    sl['TESS-p50-5cm'].ev(fGHz,theta,compensate=False)
    sl['TESS-p50-10cm'].ev(fGHz,theta,compensate=False)
    sl['TESS-p50-15cm'].ev(fGHz,theta,compensate=False)
        
    fig,ax = sl['TESS-p50-5cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='k')
    #print ax
    #fig,ax = sl['TESS-p50-10cm'].plotwrt('f',coeff='T',types=['ru'],labels=[''],color='k',linestyle='dashed',fig=fig,ax=ax)
    #fig,ax = sl['TESS-p50-15cm'].plotwrt('f',coeff='T',types=['ru'],labels=[''],color='k',linestyle='dashdot')
    plt.tight_layout()


.. image:: Glass_files/Glass_64_0.png


.. code:: python

    sl['TESS-p50-5cm'].ev(fGHz,theta,compensate=True)
    sl['TESS-p50-10cm'].ev(fGHz,theta,compensate=True)
    sl['TESS-p50-15cm'].ev(fGHz,theta,compensate=True)
    
    fig,ax = sl['TESS-p50-5cm'].plotwrt('f',coeff='T',typ=['ru'],labels=['5cm compensated',''],color='r',fig=fig,ax=ax)
    fig,ax = sl['TESS-p50-10cm'].plotwrt('f',coeff='T',typ=['ru'],labels=['10cm compensated',''],color='r',linestyle='dashed',fig=fig,ax=ax)
    fig,ax = sl['TESS-p50-15cm'].plotwrt('f',coeff='T',typ=['ru'],labels=['15cm not compensated',''],color='r',linestyle='dashdot',fig=fig,ax=ax) 
    
    fig,ax = sl['TESS-p50-5cm'].plotwrt('f',coeff='T',color='k')
    fig,ax = sl['TESS-p50-10cm'].plotwrt('f',coeff='T',color='k',linestyle='dashed',fig=fig,ax=ax)
    fig,ax = sl['TESS-p50-15cm'].plotwrt('f',coeff='T',color='k',linestyle='dashdot',fig=fig,ax=ax)

.. parsed-literal::

    /home/uguen/anaconda/lib/python2.7/site-packages/matplotlib/axes.py:4747: UserWarning: No labeled objects found. Use label='...' kwarg on individual plots.
      warnings.warn("No labeled objects found. "



.. image:: Glass_files/Glass_65_1.png


Double Glass example from litterature [1] in sub TeraHertz D-band @ 120GHz
--------------------------------------------------------------------------

.. code:: python

    sl.mat.add(name='ConcreteJc',cval=3.5,alpha_cmm1=1.9,fGHz=120,typ='THz')
    sl.mat.add(name='GlassJc',cval=2.55,alpha_cmm1=2.4,fGHz=120,typ='THz')
    sl.add('ConcreteJc',['ConcreteJc'],[0.049])
    
    theta = np.linspace(20,60,100)*np.pi/180
    sl['ConcreteJc'].ev(120,theta)
    fig,ax = sl['ConcreteJc'].plotwrt('a')
    



.. image:: Glass_files/Glass_67_0.png


.. code:: python

    plt.figure(figsize=(20,10))
    fGHz = np.linspace(110,135,50)
    sl.add('DoubleGlass',['GlassJc','AIR','GlassJc'],[0.0029,0.0102,0.0029])
    sl['DoubleGlass'].ev(fGHz,theta)
    sl['DoubleGlass'].pcolor(dB=True)


.. image:: Glass_files/Glass_68_0.png


.. code:: python

    f = plt.figure(figsize=(4,4))
    f = sl['DoubleGlass'].ev(120,theta)
    fig,ax = sl['DoubleGlass'].plotwrt('a',figsize=(10,10))
    plt.tight_layout()


.. parsed-literal::

    <matplotlib.figure.Figure at 0x7f9f03d16450>



.. image:: Glass_files/Glass_69_1.png


.. code:: python

    freq = np.linspace(110,135,50)
    sl['DoubleGlass'].ev(freq,theta)
    fig,ax = sl['DoubleGlass'].plotwrt('f',figsize=(10,10))  # @20°
    plt.tight_layout()


.. image:: Glass_files/Glass_70_0.png


References
----------

[1]. `Jacob, M. ; Kurner, T. ; Geise, R. ; Piesiewicz, R. "Reflection
ant Transmission Properties of Building Materials in D-Band for Modeling
Future mm-Wave Communication Systems" Antennas and Propagation (EuCAP),
2010 Proceedings of the Fourth European Conference
on <http://ieeexplore.ieee.org/xpl/articleDetails.jsp?tp=&arnumber=5505315&queryText%3DReflection+ant+Transmission+Properties+of+Building+Materials+in+D-Band+for+Modeling+Future+mm-Wave+Communication+Systems.QT.+Antennas+and+Propagation>`__

[2]. `R.Piesiewicz 'Terahertz characterization of building materials'
Electronics .Letters Jan 2005 Vol 41
N°18 <https://www.google.fr/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&ved=0CCwQFjAA&url=http%3A%2F%2Fwww-ece.rice.edu%2F~daniel%2Fpapers%2FnormanElecLett.pdf&ei=Tr_eUe6EG-OM0AWA0IAw&usg=AFQjCNHzt9H3RkLAtws51E9EpEgyqh-6LA&sig2=QLZlhoTJtiuHAW5Zzg_xOw&bvm=bv.48705608,d.d2k>`__

.. code:: python

    from IPython.core.display import HTML
    
    def css_styling():
        styles = open("../styles/custom.css", "r").read()
        return HTML(styles)
    css_styling()
