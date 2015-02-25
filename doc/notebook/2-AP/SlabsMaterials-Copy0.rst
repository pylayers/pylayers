
.. code:: python

    %pylab inline

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

.. code:: python

    from pylayers.antprop.slab import *
The Class ``SlabDB`` contains a dictionnary of all available Slabs. This
information is read in the file ``slabDB.ini`` of the current project
pointed by environment variable ``$BASENAME``

.. code:: python

    S = SlabDB()
A SlabDB is a dictionnary, the keys are for the current file are shown
below

.. code:: python

    S.keys()



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



Defining a new Slab and a new Material
--------------------------------------

.. code:: python

    S.mat.add(name='wall2',typ='reim',cval=2.6-0.026*1j,fGHz=4)
.. code:: python

    S.add('slab2',['wall2'],[0.15])
.. code:: python

    S.mat['wall2']



.. parsed-literal::

    {'epr': array(2.6),
     'epr2': array(-0.026),
     'epsr': (2.6-0.026j),
     'fGHz': 4,
     'index': 11,
     'mur': 1,
     'n': (1.6124717046742498-0.0080621569744854828j),
     'name': 'wall2',
     'roughness': 0,
     'sigma': 0.0057777777777777775}



.. code:: python

    S['slab2']['lmatname']



.. parsed-literal::

    ['wall2']



.. code:: python

    S['slab2']['lthick']



.. parsed-literal::

    [0.15]



.. code:: python

    fGHz= np.arange(3,5,0.01)
    theta = np.arange(0,np.pi/2,0.01)
    S['slab2'].ev(fGHz,theta)
.. code:: python

    fig = plt.figure(figsize=(10,10))
    S['slab2'].pcolor()


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_16_0.png


.. code:: python

    A=S['slab2']
As any PyLayers object there is an help function for remembering which
methods are implemented in the class.

.. code:: python

    A.help()

.. parsed-literal::

    clear: D.clear() -> None.  Remove all items from D.
    conv:  build lmat and thick
    copy: D.copy() -> a shallow copy of D
    editgui:  edit a Slab in the DB
    ev:  evaluation of the Slab
    excess_grdelay:  calculate transmission excess delay in ns
    filter:  filtering waveform
    fromkeys: dict.fromkeys(S[,v]) -> New dict with keys from S and values equal to v.
    get: D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
    has_key: D.has_key(k) -> True if D has a key k, else False
    help:  generic help
    info:  display Slab Info
    items: D.items() -> list of D's (key, value) pairs, as 2-tuples
    iteritems: D.iteritems() -> an iterator over the (key, value) items of D
    iterkeys: D.iterkeys() -> an iterator over the keys of D
    itervalues: D.itervalues() -> an iterator over the values of D
    keys: D.keys() -> list of D's keys
    loss0:  calculate loss for theta=0 at frequency (fGHz)
    losst:  Calculate loss w.r.t angle and frequency
    pcolor:  display of R & T coefficients wrt frequency an angle
    plotwrt:  plot R & T coefficients with respect to angle or frequency
    pop: D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
    popitem: D.popitem() -> (k, v), remove and return some (key, value) pair as a
    setdefault: D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D
    show:  show slab Reflection and Transmission coefficient
    tocolor:   convert slab properrties into a color
    update: D.update([E, ]**F) -> None.  Update D from dict/iterable E and F.
    values: D.values() -> list of D's values
    viewitems: D.viewitems() -> a set-like object providing a view on D's items
    viewkeys: D.viewkeys() -> a set-like object providing a view on D's keys
    viewvalues: D.viewvalues() -> an object providing a view on D's values


Information necessary to define a Slab
--------------------------------------

Each slab contains informations about its constitutive materials
electromagnetic properties.

Below an example for a simple slab, constituted with a single material
slab. The slab 'WOOD' is a layer of 4cm 'WOOD' material.

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


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_39_0.png


.. code:: python

    fGHz = np.arange(1,10,0.01)
    theta = np.arange(0,np.pi/2,0.01)
    
    S['3D_WINDOW_GLASS']['lthick']=[0.006,0.01,0.006]
    #S['3D_WINDOW_GLASS']['lmatname']=['GLASS','AIR','GLASS']
    S['3D_WINDOW_GLASS'].ev(fGHz,theta)
.. code:: python

    fig,ax = S['3D_WINDOW_GLASS'].plotwrt(var='f',coeff='T',polar='o')


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_41_0.png


.. code:: python

    fig,ax = S['WOOD'].plotwrt(var='a',coeff='R',polar='p')


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_42_0.png


plot with respect to angle

.. code:: python

    fig = plt.figure(figsize=(20,20))
    fGHz= np.array([2.4])
    S['WOOD'].ev(fGHz,theta)
    fig,ax = S['WOOD'].plotwrt(var='a',coeff='R',fig=fig)
    plt.tight_layout()



.. parsed-literal::

    <matplotlib.figure.Figure at 0x7facfc0d2c50>



.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_44_1.png


wrt to angle and frequency

.. code:: python

    plt.figure(figsize=(10,10))
    fGHz= np.arange(0.7,5.2,0.1)
    S['WOOD'].ev(fGHz,theta)
    S['WOOD'].pcolor()


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_46_0.png


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


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_47_0.png



.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_47_1.png


.. code:: python

    ## Adding new materials
.. code:: python

    theta = np.arange(0,np.pi/2,0.01)
    fGHz = np.arange(0.1,10,0.2)
    sl = SlabDB('matDB.ini','slabDB.ini')
    sl.mat.add(name='AIR2',cval=1.00000001+0j,sigma=0.00,typ='epsr')
    
    sl.add(name='AIR-5cm',lmatname=['AIR2','AIR2'],lthick=[0.05,0.05])
    sl.add(name='AIR-10cm',lmatname=['AIR2','AIR2'],lthick=[0.10,0.10])
    sl.add(name='AIR-50cm',lmatname=['AIR2','AIR2'],lthick=[0.15,0.15])
    fGHz=4
    theta = np.arange(0,np.pi/2,0.01)
    #figure(figsize=(8,8))
    # These Tessereau page 50 
    
    sl['AIR-5cm'].ev(fGHz,theta,compensate=True)
    sl['AIR-10cm'].ev(fGHz,theta,compensate=True)
    sl['AIR-50cm'].ev(fGHz,theta,compensate=True)
    
    # by default var='a' and kv = 0 
    
    fig,ax = sl['AIR-5cm'].plotwrt(color='k',labels=['5cm'])
    fig,ax = sl['AIR-10cm'].plotwrt(color='k',labels=['10cm'],linestyle='dashed',fig=fig,ax=ax)
    fig,ax = sl['AIR-50cm'].plotwrt(color='k',labels=['15cm'],linestyle='dashdot',fig=fig,ax=ax)
    plt.tight_layout()


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_49_0.png


Evaluation without phase compensation
-------------------------------------

.. code:: python

    fGHz = np.arange(2,16,0.1)
    theta = 0 
    
    sl['AIR-5cm'].ev(fGHz,theta,compensate=False)
    sl['AIR-10cm'].ev(fGHz,theta,compensate=False)
    sl['AIR-50cm'].ev(fGHz,theta,compensate=False)
    
    fig,ax = sl['AIR-5cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='r')
    #print ax
    fig,ax = sl['AIR-10cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='g',fig=fig,ax=ax)
    fig,ax = sl['AIR-50cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='b',fig=fig,ax=ax)
    sl['AIR-5cm'].ev(fGHz,theta,compensate=True)
    sl['AIR-10cm'].ev(fGHz,theta,compensate=True)
    sl['AIR-50cm'].ev(fGHz,theta,compensate=True)
    
    # by default var='a' and kv = 0 
    
    fig,ax = sl['AIR-5cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='r',linestyle='dashdot',fig=fig,ax=ax)
    fig,ax = sl['AIR-10cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='g',linestyle='dashed',fig=fig,ax=ax)
    fig,ax = sl['AIR-50cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='b',linestyle='dashdot',fig=fig,ax=ax)
    plt.tight_layout()


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_51_0.png


.. code:: python

    from pylayers.signal.bsignal import *
.. code:: python

    sl['AIR-5cm'].ev(fGHz,theta,compensate=False)
    
    S = sl['AIR-5cm']
    f=S.fGHz
    y = S.T[:,0,0,0]
    F=FUsignal(f[:,0],y)
.. code:: python

    g=F.ift(ffts=1)
.. code:: python

    g.plot(typ='v')



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7facfef992d0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7facfef999d0>]], dtype=object))




.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_55_1.png


.. code:: python

    sl['AIR-5cm'].ev(fGHz,theta,compensate=True)
    sl['AIR-10cm'].ev(fGHz,theta,compensate=True)
    sl['AIR-50cm'].ev(fGHz,theta,compensate=True)
    
    fig,ax = sl['AIR-5cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='k')
    #print ax
    fig,ax = sl['AIR-10cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='k',linestyle='dashed',fig=fig,ax=ax)
    fig,ax = sl['AIR-50cm'].plotwrt('f',coeff='T',typ=['ru'],labels=[''],color='k',linestyle='dashdot',fig=fig,ax=ax)
    plt.tight_layout()



.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_56_0.png


.. code:: python

    sl.mat.add(name='ConcreteJc',cval=3.5,alpha_cmm1=1.9,fGHz=120,typ='THz')
    sl.mat.add(name='GlassJc',cval=2.55,alpha_cmm1=2.4,fGHz=120,typ='THz')
    sl.add('ConcreteJc',['ConcreteJc'],[0.049])
    
    theta = np.linspace(20,60,100)*np.pi/180
    sl['ConcreteJc'].ev(120,theta)
    fig,ax = sl['ConcreteJc'].plotwrt('a')
    



.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_57_0.png


.. code:: python

    plt.figure(figsize=(20,10))
    fGHz = np.linspace(110,135,50)
    sl.add('DoubleGlass',['GlassJc','AIR','GlassJc'],[0.0029,0.0102,0.0029])
    sl['DoubleGlass'].ev(fGHz,theta)
    sl['DoubleGlass'].pcolor(dB=True)


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_58_0.png


.. code:: python

    f = plt.figure(figsize=(4,4))
    f = sl['DoubleGlass'].ev(120,theta)
    fig,ax = sl['DoubleGlass'].plotwrt('a',figsize=(10,10))
    plt.tight_layout()


.. parsed-literal::

    <matplotlib.figure.Figure at 0x7fad016eda10>



.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_59_1.png


.. code:: python

    freq = np.linspace(110,135,50)
    sl['DoubleGlass'].ev(freq,theta)
    fig,ax = sl['DoubleGlass'].plotwrt('f',figsize=(10,10))  # @20°
    plt.tight_layout()


.. image:: SlabsMaterials-Copy0_files/SlabsMaterials-Copy0_60_0.png


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
