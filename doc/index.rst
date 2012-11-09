
PyLayers : Python radio propagation and  localization simulator 
===============================================================

Pylayers is a simulator designed to model mobile pedestrian radio propagation.
Its purpose is to help designing localization and tracking algorithms for advanced 
communications systems. Hopefully, also useful for
higher education on topics related to wireless communications. 


Example 
=======

.. ipython::
    

    In [1]: import pylayers.simul.simulem as sem

    In [1]: S = sem.Simul()


Features
========

* Layout edition
* Extensive graph description of layout  
* Full space vector antenna patterns using a sparse spherical harmonics representation  
* Multi layers transmission coefficients 
* Motley Keenan model 
* Ray signature determination and analysis 
* Embedded C coded UWB ray tracing engine
* Handling of Impulse Radio Ultra Wideband Waveforms 
* Mobile user mobility based on E 
* Set Membership positoning algorithm toolbox
* Heterogeneous positioning toolbox 
* Various techniques exploiting heterogeneous Radio observables



Context
=======

PyLayers is beeing developped for research into radio localization exploiting impulse radio
signals as being defined in the IEEE 802.15.4a and 802.15.6 standards.
Pylayers could also evolved to help with studies related to channel characterization,
Wireless Sensor Network, MANET, radio cooperation and relaying techniques, 
physical layer security, WBAN applications, advanced home automation,
connected sensito-cognitive environments ...

Documentation
=============

.. toctree::
    :maxdepth: 1 

    manual/index.rst 

    modules/index.rst
    
    bibliography.rst



Third party packages 
====================

Pylayers is written in Python and tributes highly the following Python packages

* numpy ( for algebra in multi dimensions )  
* scipy ( for scientific libraries )
* simpy ( for discrete events simulation ) 
* matplotlib ( for the graphical outputs ) 
* networkx ( for graph description and algorithms ) 
* shapely ( for planar geometry )   
* scikit-learns( for cutting edge machine learning tools ) 
* cvxopt ( for convex optimization ) 
* spherepack (as spherical harmonics library)  
* ipython ( to rule them all ... ) 


In the near future probably for acceleration of critical part of the code :  

* Cython 
* numba 

Pylayers relies also on the good old `geomview` tool for fast and simple interaction with 3D entities. 


All those dependencies make the installation of pylayers somehow tricky. 

`Pylayers` is modular and designed to be easily to extend, 
allowing it to evolve over the time as a tool for research and development in a wide range 
of communications and localization topics.

We warmly encourage all users to contribute 
new algorithms, models and other improvements back to the project.

**Just fork it !** on your github account.


Developpers
===========

Pylayers is driven by professor Bernard Uguen at `University of Rennes 1
<http://www.univ-rennes1.fr/>`_, `IETR laboratory <http://www.ietr.fr/>`_ and 
`ESIR school of engineering <http://esir.univ-rennes1.fr/>`_

Pylayers is currently developed at IETR  by Bernard Uguen, Nicolas Amiot, Mohamed
Laaraiedh and Meriem Mhedbhi, with the technical support of all the members from the
Research Team of the Propagation and Localization team of IETR lab.    

Early contributors : Friedman Tchoffo Talom, Roxana Burghelea, Yu Lei, Taguhi
Chaluyman, Stéphane Avrillon and Eric Plouhinec.

Acknoledgments 
==============

This work has been supported by the Bretagne Region council (Project LOCUS), by the french ANR
project AUBADE and CORMORAN and by the European projects `FP7 UCELLS <http://www.ist-ucells.org/>`_, 
FP7 WHERE1, `FP7 WHERE2 <http://www.ict-where.eu/>`_.

.. image:: _static/logoUR1.jpg 
    :scale: 18%

.. image:: _static/logoIETR.jpg 
    :scale: 22%

.. image:: _static/logoESIR.png 
    :scale: 20%

.. image:: _static/bretagnegd.jpg 
    :scale: 20%

.. image::  _static/ucells.png 
    :scale: 40%

.. image:: _static/where1logo.jpg 
    :scale: 30%

.. image:: _static/WHERE2_Logo.jpg 
    :scale: 20%

.. image:: _static/logo_CORMORAN.png 
    :scale: 20%

.. image:: _static/IR.png 
    :scale: 30%


Download and Installation
-------------------------

The current version is tagged 0.1, released on 9 November 2012. Download the
release in your preferred format on github.

.. code-block:: bash
    
    $ git clone https://github.com/buguen/pylayers.git


Mailing List
------------


There is a "http://groups.google.com/group/pylayers-users" mailing list for
users of Pylayers on which you can get help with using or extending the
software. New releases will also be announced on the mailing list. 


License
-------

Copyright © 2012 University of Rennes 1

Pylayers is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Pylayers is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

