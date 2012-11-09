
PyLayers : Python radio propagation and  localization simulator 
===============================================================

Pylayers is a simulator designed to model mobile pedestrian radio propagation.
The purpose is to provide a tool helping the design of localization and tracking
algorithms for advanced communications systems.


Example 
=======

.. ipython::
    
    In [1]: import pylayers.simul.simulnet as snet

    In [1]: import pylayers.simul.simulem as sem

    In [1]: Sn = snet.Simul()

    In [1]: Se = sem.Simul()


Features
========

* Indoor layout description and edition
* Extensive graph description of indoor layout  
* Multi layers transmission coefficients 
* Account fot the full vector antenna pattern  
* Motley Keenan model 
* Ray Signature determination and analysis 
* Embedded C coded modular ray tracing engine
* Handling of Impulse Radio Ultra Wideband Waveforms 
* Mobile user mobility 
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

The simulator is written mostly in Python and depends on the following packages

* numpy
* scipy
* simpy
* matplotlib
* networkx   
* shapely   
* scikit-Learn
* cvxopt
  
`PYLAYERS` is designed to be easy to extend with new functionalities, 
allowing it to evolve over the time as a common
tool for research and development in a wide range of communications and
localization applications.

We warmly encourage all users to contribute 
new algorithms, models and other improvements back to the project.

Just fork the github repository on github.

Pylayers was written by Bernard Uguen, Nicolas Amiot, Mohamed Laaraiedh, Yu Lei, Roxana Burghelea
, Friedman Tchoffo Talom at the IETR, University of Rennes 1.

This work has been supported by the Bretagne Region council (Project LOCUS), by the french ANR
project AUBADE and CORMORAN and by the European projects FP7 UCELLS, FP7 WHERE1, FP7 WHERE2.


Documentation
-------------

.. toctree::
    :maxdepth: 1 

    manual/index.rst 

    modules/index.rst
    
    bibliography.rst


Download and Installation
-------------------------

The current version is tagged 0.1, released on 9 November 2012. Download the release in your preferred format:

.. topic::
    
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

