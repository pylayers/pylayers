
PyLayers : Python radio propagation and  localization simulator 
===============================================================

Pylayers is a simulator designed to model mobile pedestrian radio propagation.
Its purpose is to help the proper design of localization and tracking
algorithms in advanced communications systems. Hopefully also useful for
higher education on topics related to wireless communications. 


Example 
=======

.. ipython::
    

    In [1]: import pylayers.simul.simulem as sem

    In [1]: S = sem.Simul()


Features
========

* Indoor layout description and edition
* Extensive graph description of indoor layout  
* Account for the full vector antenna patterns using a sparse spherical harmonics representation  
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

PyLayers has been developped to help with research into impulse radio
propagation studies and  algorithms for
localization, but it may be also useful to help with the development of
characterization, Wireless Sensor Network, MANET, radio cooperation techniques, physical layer
security ...


Third party packages 
====================

Pylayers is written in Python and tributes highly the following Python packages

* numpy ( for the algebra in higher dimensions and just for the fun of using `einsum` function )  
* scipy ( for the huge amount of science it contains )
* simpy ( for elegant discrete events simulation ) 
* matplotlib ( for the Art ) 
* networkx ( for the elegance of graph theory ) 
* shapely ( for the planar geometry )   
* scikit-Learn ( for the cutting edge machine learning tools ) 
* cvxopt ( for convex optimization ) 
* ipython ( to rule them all ... ) 


In the future probably also :  
* Cython 
* numba 

Pylayers relies also on the good old `geomview` tool for fast and simple interaction with 3D entities. 


All those dependencies make the installation of pylayers somehow tricky. This
is a good 
`PYLAYERS` is modular and designed to be easily to extend, 
allowing it to evolve over the time as a tool for research and development in a wide range 
of communications and localization topics.

We warmly encourage all users to contribute 
new algorithms, models and other improvements back to the project.

**Just fork it !** on your github account.


Developpers
===========

Pylayers is driven by professor Bernard Uguen the `University of Rennes 1 <http://www.univ-rennes1.fr/>`_
(`IETR laboratory <http://www.ietr.fr/>`_ and `ESIR school of engineering <http://esir.univ-rennes1.fr/>`_)
The code is currently developed at IETR  by Bernard Uguen, Nicolas Amiot, Mohamed
Laaraiedh and Meriem Mhedbhi, with the technical support of all the members from the
Research Team of the Propagation and Localization team of IETR lab.    

Early contributors : Yu Lei, Roxana Burghelea, Friedman Tchoffo Talom,  Stéphane Avrillon and Eric Plouhinec.

This work has been supported by the Bretagne Region council (Project LOCUS), by the french ANR
project AUBADE and CORMORAN and by the European projects `FP7 UCELLS <http://www.ist-ucells.org/>`_, 
FP7 WHERE1, `FP7 WHERE2 <http://www.ict-where.eu/>`_.

.. image:: _static/logoIETR.jpg 
    :scale: 20%

.. image:: _static/bretagne.png 
    :scale: 20%

.. image::  _static/ucells.png 
    :scale: 20%

.. image:: Where1.png 
    :scale: 20%

.. image:: Where2.png 
    :scale: 20%

.. image:: Cormoran.png 
    :scale: 20%

Documentation
-------------

.. toctree::
    :maxdepth: 1 

    manual/index.rst 

    modules/index.rst
    
    bibliography.rst


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

