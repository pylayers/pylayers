Simulation File        
---------------

This example load a simulation file and shows the layout with the surimposed
grid of tx and rx points which define the links of interest for the simulation 

.. ipython::
    
    In [1]: from pylayers.simul.simulem import *
    In [1]: from pylayers.gis.layout import *
    In [1]: from numpy import *
    In [1]: import matplotlib.pylab as plt  

    In [1]: S = Simul('default.ini')

    @savefig DefaultLayout.png width=8in 
    In [1]: S.L.showGs()


Getting started 
---------------

First of all, it is necessary to fill an **.ini** file which gathers
informations required for starting a simulation. If no file 

Simulation `.ini`  file
~~~~~~~~~~~~~~~~~~~~~~~

Below is presented an example of a simulation file `default.ini`::


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
--------------

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


