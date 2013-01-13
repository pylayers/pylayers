Setting a first simulation file        
--------------------------------

An example simulation file is available in the directory which contains
`pylayers/data/ini`

This **.ini** file gathers informations required for starting a simulation. If no file 
is present a default file is created. 

Simulation `default.ini`  file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example of a simulation file `default.ini`::


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
                list of calculated output files



Output Section 
--------------

The output section is used to keep track of calculated links. The key 
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

This example loads a simulation file and shows the layout with the surimposed
tx and rx grid points which defines the links of interest for the simulation 

.. ipython::
    
    In [1]: from pylayers.simul.simulem import *
    In [1]: from pylayers.gis.layout import *
    In [1]: from numpy import *
    In [1]: import matplotlib.pylab as plt  

    In [1]: S = Simul('def.ini')

    @savefig DefaultSimul.png width=8in 
    In [1]: S.show()


.. ipython::

    @savefig DefaultLayout.png width=8in 
    In [1]: S.L.showGs()

