
.. code:: python

    from IPython.core.display import HTML
    from IPython.display import FileLink
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

    import ConfigParser
    import pylayers.util.pyutil as pyu

Network Simulation Configuration
================================

PyLayers is designed to provide indoor radio channel simulation for
mobile agents.

The goal is to address mobility in indoor environment heterogeneous
network, with human being carriers of a mobile User Equipement (UE)
which possibly embeds several Radio Acess Technology (RAT).

Several humans can be created and their motion in the environement
should be as realistic as possible, because for many applications it
turns out that many parameters of interest are stongly dependent of the
dynamic topology of the mobile network.

In the following the configuration files for proceeding with those high
level ``PyLayers`` simulation are described.

The configuration file is ``simulnet.ini``

Simulnet.ini
------------

This file is assumed to be located in ``$BASENAME/ini``. As the format
of this file is not stable yet. Refer to the comment in the file below
for obtaining the inforamtion about the format.

.. code:: python

    !cat $BASENAME/ini/simulnet.ini

.. parsed-literal::

    [Mysql]
    host = localhost
    user = root
    passwd = sqlsql
    dbname = test
    dumpdb =True
    
    [Save]
    ; deprecated save option
    save=[]
    ;save=['csv','mysql','matlab','pyray','txt','ini']
    ; save format using Save class..to be deprecatred soon
    savep=True
    ; pandas save format . only record mechanical. To be upgraded soon
    savepd=True
    
    [Layout]
    filename = TA-Office.ini
    
    x_offset  = 30
    y_offset = 2
    
    the_world_width	 = 65
    the_world_height = 20
    the_world_scale	 = 20 
    
    [Mechanics]
    ; update time for agent movement
    mecanic_update_time = 0.2
    ; select how agnt choose destiantion
    ;'random' ; file
    choose_destination = 'random'
    ; experimental show for debug purpose
    pdshow=False
    
    [Network]
    ; simulate the network
    network=True
    ; refresh TOA regulary 'synchro 'or with distance 'autionomous'
    Communication_mode='autonomous'
    ; update time for refreshing network
    network_update_time = 0.1
    ; show nodes moving & radio link
    show = False
    ; show in ipython notebook
    ipython_nb_show = False
    ; show signature ( not fully functionnal)
    show_sg = False
    ; show 2 tables : mecanic & network
    show_table = False
    ; show the same information but in terminal
    dispinfo = False
    
    [Localization]
    ; perform localization
    localization = False
    ; time to refresh localization
    localization_update_time = 1.0
    ; list of used methods method = ['alg','geo']
    method = ['geo']
    
    
    
    [Simulation]
    ; simulation filename for savepd
    filename = 'simulnet'
    ; Simulation duration
    duration = 120.0
    ; speed ratio ag
    speedratio = 50.
    ; time for refreshing tk plot ( obsolete)
    show_interval = 0.5
    ; show scene using tk renderer ( obsolete)
    showtk   = False
    ; choose seed for random mobiliity
    seed = 1
    ; verbose output
    verbose = False
    


.. code:: python

    Cp = ConfigParser.ConfigParser()
    Cp.read(pyu.getlong('simulnet.ini','ini'))



.. parsed-literal::

    ['/home/uguen/Bureau/P1/ini/simulnet.ini']



Current version of ``Simulnet.ini`` contains the following sections

.. code:: python

    Cp.sections()



.. parsed-literal::

    ['Mysql',
     'Save',
     'Layout',
     'Mechanics',
     'Network',
     'Localization',
     'Simulation']



Save section
~~~~~~~~~~~~

The save section handles the output files of the simulation.

.. code:: python

    dict(Cp.items('Save'))



.. parsed-literal::

    {'save': '[]', 'savep': 'True', 'savepd': 'True'}



The ``savep`` boolean enable/disable saving of the simulation.

.. code:: python

    dict(Cp.items('Save'))['savep']



.. parsed-literal::

    'True'



The log file which contains all traces from the simulated dynamics are
in ``$BASENAME/netsave``

.. code:: python

    !ls $BASENAME/netsave/*

.. parsed-literal::

    /home/uguen/Bureau/P1/netsave/save  /home/uguen/Bureau/P1/netsave/save.mat  /home/uguen/Bureau/P1/netsave/simulnet_TA-Office.h5  /home/uguen/Bureau/P1/netsave/traj_nicta.h5


Layout section
~~~~~~~~~~~~~~

This section specifies the layout parameter and spatial dimension of the
simulation

.. code:: python

    dict(Cp.items('Layout'))



.. parsed-literal::

    {'filename': 'TA-Office.ini',
     'the_world_height': '20',
     'the_world_scale': '20',
     'the_world_width': '65',
     'x_offset': '30',
     'y_offset': '2'}



Choose the used Layout for simulation

.. code:: python

    dict(Cp.items('Layout'))['filename']



.. parsed-literal::

    'TA-Office.ini'



Setup an offset for defining the coordinate system origin

.. code:: python

    print dict(Cp.items('Layout'))['x_offset']
    print dict(Cp.items('Layout'))['y_offset']

.. parsed-literal::

    30
    2


Network section
~~~~~~~~~~~~~~~

.. code:: python

    dict(Cp.items('Network'))



.. parsed-literal::

    {'communication_mode': "'autonomous'",
     'dispinfo': 'False',
     'ipython_nb_show': 'False',
     'network': 'True',
     'network_update_time': '0.1',
     'show': 'False',
     'show_sg': 'False',
     'show_table': 'False'}



Setup communication mode between node:

-  ``"autonomous"`` : the data exchange between nodes is driven by the
   localization layer. If more information is required to estimate the
   position then a communication request is sent to the communication
   state
-  ``"synchro"`` : the data exchange between nodes is periodic. LDPs are
   periodically refreshed at the ``network_update_time``

.. code:: python

    dict(Cp.items('Network'))['communication_mode']



.. parsed-literal::

    "'autonomous'"



Time step for the refresh network information

.. code:: python

    dict(Cp.items('Network'))['network_update_time']



.. parsed-literal::

    '0.1'



Vizualization of the simulation using matplotlib

.. code:: python

    dict(Cp.items('Network'))['show']



.. parsed-literal::

    'False'



Vizualization of a table summing up the data exchange of the nodes

.. code:: python

    dict(Cp.items('Network'))['show_table']



.. parsed-literal::

    'False'



Vizualization of the simulation inside ipython notebook

.. code:: python

    dict(Cp.items('Network'))['ipython_nb_show']



.. parsed-literal::

    'False'



Mechanics
---------

This section specifies agents dynamic during simulation

.. code:: python

    dict(Cp.items('Mechanics'))



.. parsed-literal::

    {'choose_destination': "'random'",
     'mecanic_update_time': '0.2',
     'pdshow': 'False'}



Setup how agent choose their target:

-  ``"random"``: the agnet move into the layout randomly
-  ``"file"`` : the agent follow the sequence specified in
   ``<project_dir>/nodes_destination.ini``

.. code:: python

    dict(Cp.items('Mechanics'))['choose_destination']



.. parsed-literal::

    "'random'"



Time step for refreshing the mechanical layer (ground truth position)

.. code:: python

    dict(Cp.items('Mechanics'))['mecanic_update_time']



.. parsed-literal::

    '0.2'



Localization section
~~~~~~~~~~~~~~~~~~~~

Setup Localization algorithms

.. code:: python

    dict(Cp.items('Localization'))



.. parsed-literal::

    {'localization': 'False',
     'localization_update_time': '1.0',
     'method': "['geo']"}



enable/disable localizaiton of the agents

.. code:: python

    dict(Cp.items('Localization'))['localization']



.. parsed-literal::

    'False'



Select localization methods :

-  Algebraic : hétérogeneous localization algorithm
-  Geometric : RGPA

.. code:: python

    dict(Cp.items('Localization'))['method']



.. parsed-literal::

    "['geo']"



Time step for localization update

.. code:: python

    dict(Cp.items('Localization'))['localization_update_time']



.. parsed-literal::

    '1.0'



Simulation section
~~~~~~~~~~~~~~~~~~

.. code:: python

    dict(Cp.items('Simulation'))




.. parsed-literal::

    {'duration': '120.0',
     'filename': "'simulnet'",
     'seed': '1',
     'show_interval': '0.5',
     'showtk': 'False',
     'speedratio': '50.',
     'verbose': 'False'}



Setup simulation duration in second

.. code:: python

    dict(Cp.items('Simulation'))['duration']




.. parsed-literal::

    '120.0'



Setup random seed for simulation

.. code:: python

    dict(Cp.items('Simulation'))['seed']




.. parsed-literal::

    '1'



Display messages during simulation

.. code:: python

    dict(Cp.items('Simulation'))['verbose']




.. parsed-literal::

    'False'



See Also

.. code:: python

    FileLink('../4-MOB/Mobility.ipynb')



.. raw:: html

    <a href='../4-MOB/Mobility.ipynb' target='_blank'>../4-MOB/Mobility.ipynb</a><br>


