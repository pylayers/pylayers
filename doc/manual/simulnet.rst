How to goes from a simulnet to a Channel Impulse Response ? 
===========================================================

Generate a mobile simulation with simulnet
-------------------------------------------

You need first to have an `agent.ini` and `simulnet.ini` file in the  `/ini` directory
of your current project. 

agent.ini 
~~~~~~~~~

Below is an exemple of agent.ini file 

.. literal::

    [used_agent]
    list=['A1','A2','BS1','BS2']
    ;list=['A1','A2']

    [A1]
    name = John
    ID = 1
    type = ag
    roomId = 0
    pos=[]
    RAT = ['rat1']
    epwr = [0]
    refreshRSS = 0.3
    refreshTOA = 0.8

    condition= "distance:'<5';node:'all';rat:'all';message:['d','pe'];and"
               "distance:'>10';node:'all';rat:'rat2';message:['pe'];"

    [A2]
    name = Steve
    ID = 2
    type = ag
    roomId = 15
    pos=[]
    RAT = ['rat1']	
    epwr = [0]
    refreshRSS = 0.4
    refreshTOA = 0.7

    condition= "distance:'<5';node:'all';rat:'rat1';message:['d','pe'];and"
               "distance:'>10';node:'all';rat:'rat2';message:['pe'];"

    [A3]
    name = Peter
    ID = 3
    type = ag
    roomId = 2
    pos=[]
    RAT = ['rat1','rat2']	
    epwr = [0,0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    condition= "distance:'<5';node:'all';rat:'rat1';message:['d','pe'];and"
               "distance:'>10';node:'all';rat:'rat2';message:['pe'];"


    [A4]
    name = Mike
    ID = 4
    type = ag
    roomId = 8
    pos=[]
    RAT = ['wifi','bt','lte']	
    epwr = [0,0,0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    [A5]
    name = Micheal
    ID = 5
    type = ag
    roomId = 15
    pos=[]
    RAT = ['wifi','bt','lte']	
    epwr = [0,0,0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    [BS1]
    name = AP1
    ID = 6
    type = ap
    roomId = -1
    ;pos= [-5,5]
    pos= [0.5,2.]
    RAT = ['rat1']	
    epwr = [0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    condition= "distance:'<5';node:'all';rat:'rat1';message:['d','pe'];and"


    [BS2]
    name = BS2
    ID = 7
    type = ap
    roomId = -1
    ;pos= [70,5]
    pos= [0.7,14]
    RAT = ['rat1']	
    epwr = [0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    condition= "distance:'<5';node:'all';rat:'rat1';message:['d','pe'];and"

    [BS3]
    name = BS3
    ID = 8
    type = ap
    roomId = -1
    ;pos= [0,20]
    pos= [39.,13.]
    RAT = ['rat1']
    epwr = [0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    condition= "distance:'<5';node:'all';rat:'rat1';message:['d','pe'];and"

    [BS4]
    name = BS4
    ID = 9
    type = ap
    roomId = -1
    pos= [70,20]
    RAT = ['rat1']	
    epwr = [0]
    refreshRSS = 0.5
    refreshTOA = 0.5

    condition= "distance:'<5';node:'all';rat:'rat1';message:['d','pe'];and"


configure simulnet.ini file

.. literal::

    [Mysql]
    host = localhost
    user = root
    passwd = sqlsql
    dbname = test
    dumpdb =True

    [Save]
    save=['pyray','txt']
    ;save=['csv','mysql','matlab','pyray','txt','ini']


    [Layout]
    filename = TA-Office.str

    x_offset  = 30
    y_offset = 2

    the_world_width	 = 65
    the_world_height = 20
    the_world_scale	 = 20 

    [Mecanic]
    ; update time for agent movement
    mecanic_update_time = 0.1

    [Network]
    ; update time for refreshing network
    network_update_time = 0.1
    ; show nodes moving & radio link
    show = False
    ; show signature ( not fully functionnal)
    show_sg = False
    ; show 2 tables : mecanic & network
    show_table = False
    ; show the same information but in terminal
    dispinfo = False

    [Localization]
    ; not implemented yet
    ; perform localization
    localization = True
    ; time to refresh localization
    localization_update_time = 0.25

    [Simulation]
    ; Simulation duration
    duration = 1.0
    ; time for refreshing tk plot ( obsolete)
    show_interval = 0.5
    ; show scene using tk renderer ( obsolete)
    showtk   = False



`[Save]` section 
-----------------

+ `pyray` : use this option to generate outputs of the simulation compliant with pulsray  
+ `txt` : use this option to generate outputs of the simulation compliant with the WHERE2 DB

`[Layout]` section :
-------------------

+ `filename` option give the name ofd the file .str used for simulnet simulation

`[Network]` section 
-------------------

+ The network_update_time option give the sample rate of the output files.
( ! For a correct output GIVE THE SAME VALUE TO
mecanic_update_time option in the [Mecanic] section)
+ The show option allow to display with matplotlib the simulation trace.( True
or False)
In section [Localization] :
+ localization : boolean True/false for each mobile node (agent) compute their
position
+ localization_update_time : refresh time for localization
In Section [Simulation]:
+ The duration option set up the simulation duration

configure agent.ini file

Before running the simulation you can select the involved agents into the
simuation.
Go to ProjectDirectory/ini/agent.ini
In Section [used_agent]:
+ append the list option with the list of available agent in the following of
the file . Ex: list=['A1','A2','BS1','BS2']
In Section [agentidentifier]:

+ name name of node
+ ID identifier of the node
+ type option : select ag for agent ( mobile node) or ap for acces point
( static node)
+ roomId room where the node start the simulation
+ pos if roomID=-1 : position [x,y] where the node start the simulation
+ RAT option : select a list of rat name ( for now just keep only 'rat1') for
+ epwr list of emmitted power for each rat
+ refreshRSS refresh frequency of RSS measurement for the node
+ refreshTOA refresh frequency of TOA measurement for the node
+ condition not used for now


Running the simulation
-----------------------

.. python::

    >>> import pylayers.simul.simulnet as snet 
    >>> S = snet.Simul()
    >>> S.runsimul()

Illustration 1: terminal output after S.runsimul()

+ would you like to erase previous txt files ?
+ answer y to replace the previous .txt files ( compliant with the W2
database) by the ones from the current simulation


Once the simulation is finished, generated files can be find into

ProjectDirectory/netsave

Illustration 2: contents of ProjectDirectory/netsave


on the illusration 2

+ 1 and 2 are 2 mobile nodes. 1.ini and 2.ini contains trajectories compliant with pyray
+ 6 and 7 are 2 anchors . 6.ini and 7.ini contains a single position compliant with pyray
+ pyray.ini contains some information from the simulnetsimulation
+ all .txt are generated for the W2 database

Illustration 3: exemple of contents of 1.ini file


Compute and exploit pyray with the given simulnet files
-------------------------------------------------------

computation

.. python::

    >>> from pylayers.simul.exploit import *
    >>> E=Exploit()
    >>> E.compute()


E.compute runs the pyray simulation for
+ all anchors nodes (ap) to all mobile nodes (ag)
+ all mobile nodes (ag) a time 't' all mobile nodes (ag) at the same time 't'

This operation can take a while.....

Results
-------

The simulation generates .mat files in directory `ProjectDirectory/output/nodeid`

Illustration 4: Contents of ProjectDirectory/output/1 and 6

file name is build as follow : `defaultcir-tx_node_id-rx_node_id-pposition_id`

The position_id is related to the position number into
`ProjectDirectory/netsave/nodeid.ini`

Exploitation
------------

either

.. python::

    >>> E.pltcir( nodeid1 , nodeid2, position_number)
or

.. python::

    >>> E.pltciri( nodeid1 , nodeid2)

 With pltciri you can interract with the plot to display the desired CIR.

Invocate pltciri for 2 nodes for instance 1 and 6:

.. python::

    >>> E.pltciri(1,6)

The Layout and the nodes tracjectories are diplayed:

On the figure press 't' and then click a point of the node trajectory ( for instance the
red dot), to select a transmitter Tx

On the figure press 'x' and then click a point of the node trajectory ( for instance the
one of the blue) to select a receiver Rx
Press Enter to show the associated Channel impulse response (CIR)

+ Selected dot are marked with a black cross
