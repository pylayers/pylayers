module simulem
==============

The module simulem is in charge of calculating the received UWB waveform for a
given radio link. This section is built from the file test_simulem.py in the
tests sub-directory of antprop directory. 


.. warning:
    This demo works only on LInux
    platform. It relies on exec files (launching, tracing,tratotud,evalfield) in the `bin` 
    directory of the project. It is foreseen to replace soon those dependencies by
    pure python equivalent. 

Simulation File        
---------------

This example load a simulation file and shows the layout with the surimposed
grid of tx and rx points which define the links of interest for the simulation 



Getting started 
---------------

First of all, it is necessary to fill an **.ini** file which gathers
informations required for starting a simulation. If no file exists one 
is default file is created. 

structure of the simulation description file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is presented an example a `default.ini` simulation file ::


    [files]
    mat = matDB.ini
    tx = radiotx.ini
    slab = slabDB.ini
    txant = defant.vsh3
    rx = radiorx.ini
    patra = def.patra
    conf = project.conf
    palch = def.palch
    struc = Lstruc.str
    rxant = defant.vsh3

    [waveform]
    tw = 30
    band = 0.499
    fc = 4.493
    thresh = 3
    fe = 50
    type = generic

    [frequency]
    fghzmin = 2.0
    fghzmax = 11.0
    nf = 181

    [tud]
    purc = 100
    num = -1
    nrmax = 500

    [output]
    1 = default1.ini

    

This file is composed of independant sections which are respectively ::

        [files]
                This section contains the short name of the required input file 
                for high level commands 
        [launching]
                various parameters for the launching phase 
        [tracing]
                various parameters for the tracing phase 
        [waveform]
                parameters defining the applied waveform 
        [frequency]
                electromagnetic frequency range 
        [tud]
                ray filtering parameters
        [output]
                already calculated output files



Output Section 
~~~~~~~~~~~~~~

The output section is used to keep track of already calculated links. The key 
is an integer which correspond to a radionode index and the corresponding
associated value is a file which is stored in the `output` directory of the
project. 

Below is an example of the content of an output `.ini` file ::


    [rang]
    1 = defstr_slabDB_def_radiotx_1_def_radiorx_1_0_500.rang

    [trace]
    1 = defstr_slabDB_def_radiotx_1_def_radiorx_1.tra

    [launch]
    1 = defstr_slabDB_def_radiotx_1.lch

    [tang]
    1 = defstr_slabDB_def_radiotx_1_def_radiorx_1_0_500.tang

    [tauk]
    1 = defstr_slabDB_def_radiotx_1_def_radiorx_1_0_500.tauk

    [field]
    1 = defstr_slabDB_def_radiotx_1_def_radiorx_1_0_500.field

    [tud]
    1 = defstr_slabDB_def_radiotx_1_def_radiorx_1_0_500.tud

    [cir]
    1 = cir-tx001-rx001


Example
-------

Modules are loded 

.. ipython::
    
    In [1]: from pylayers.simul.simulem import *
    In [1]: from pylayers.gis.layout import *
    In [1]: from numpy import *
    In [1]: import matplotlib.pylab as plt  

A simulation is loaded from file 'default.ini' 

.. ipython::

    In [1]: S = Simul('defstr.ini')

    @savefig DefaultLayout.png width=8in 
    In [1]: S.L.showGs()
A layout file is associated to simulation and loaded and displayed.

.. code:: python 

    filestr = 'defstr'
    S.layout(filestr+'.str','matDB.ini','slabDB.ini')
    S.L.showGs()
    plt.show()


.. code::  python     
    setting transmitter
    S.tx = RadioNode(typ='tx')
    S.tx.point([1.2,1,1.4])

    # setting receiver

    S.rx = RadioNode(typ='rx')
    S.rx.point([8,-1.2,1.5])
    S.rx.point([8,-1.21,1.5])
    S.rx.point([8,-1.22,1.5])
    S.rx.point([8,-1.23,1.5])
    S.rx.point([8,-1.24,1.5])

    # saving simulation 

    S.save()

    # print launching parameters
    S.palch.info()

    # ang Tx : angular step from Tx
    S.palch.angTx  = 1

    # ISB ang Incident Shadow Boundary angle (degree) 
    S.palch.ISBang = 90  

    # ray elimination Threshold 
    S.palch.ethreshold = 0.001

    # maximum depth
    S.palch.maxdeep  = 10

    # typealgo = 0 (include diffraction) 1 (no diffraction)
    S.palch.typalgo = 1
    title = str(S.palch.angTx) + '-' +\
            str(S.palch.ISBang) + '-' +\
            str(S.palch.ethreshold) + '-' + \
            str(S.palch.maxdeep) + '-' + \
            str(S.palch.typalgo)

    S.palch.save()
    S.pafreq.fghzmin=2
    S.pafreq.fghzmax=11
    S.pafreq.nf=181
    S.pafreq.save()
    # showing the simulation 
    print "Launching "
    print "-----------------"
    S.launching(1)

    # retrieve the launching tree

    L1 = S.getlaunch(1)

    # display the launching tree for different depths

    fig = plt.figure(figsize=(10,10))
    plt.title('launching parameters '+title+' '+filestr )
    plt.axis('off')
    N = S.palch.maxdeep
    M = N/2
    #
    #for k in range(N):
    #    ax = fig.add_subplot(M,2,k+1)
    #    fig,ax = L1.show(S.L,k+1,f=fig)

    #    fig.savefig(pylayersdir+'/doc/auto_examples/simul/'+filestr+'-launching.png')    
    print "Tracing "
    print "-----------------"
    print "purc :",S.config.get('tud','purc')
    fig = plt.figure()
    S.tracing(1,1)
    gr = GrRay3D()
    gr.load(S.dtra[1][1],S.L)
    #f,a = S.L.showGs(fig=fig)
    #plt.axis('on')
    #gr.show(fig=f,ax=a,rayset=np.arange(100))
    print "Tratotud "
    print "-----------------"
    print "purc :",S.config.get('tud','purc')
    S.tratotud(1,1)
    gt = GrRayTud()
    # loading rays in tud format 
    #gt.load(S.dtud[1][1],S.dtang[1][1],S.drang[1][1],S.L.sl)
    #print "Evalfield "
    #print "-----------------"
    S.field(1,1)
    S.cir(1,1)
    f = plt.figure()
    S.pltcir(1,1,fig=f)
    plt.show()
