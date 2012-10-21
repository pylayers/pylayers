PYLAYERS Overview
=================

Pylayers is a simulator designed to model mobile pedestrian radio propagation
in Indoor environments.


The simulator is able to welcome various kind of radio propagation solvers of 
unequal complexity. The current version can generate UWB channel
impulse responses along trajectory thanks to external C functions which have
been developed few years ago during the french ANR project AUBADE.


Taking advantage of the project FP7 UCELLS, FP7 WHERE1 and FP7 WHERE2 the tool
has gradually been dedicated to localization applications. Pylayers contains
various functions which allows to estimate the position based a given number
of Location Dependent Parameters (LDPs)

Features included in the current development version are:

+ Indoor layout description and edition
+ Extensive graph description of indoor layout  
+ Multi layers transmission coefficients 
+ Account fot the full vector antenna pattern  
+ Motley Keenan model 
+ Ray Signature determination and analysis 
+ Embedded C coded modular ray tracing engine
+ Handling of Impulse Radio Ultra Wideband Waveforms 
+ Mobile user mobility 
+ Set Membership positoning algorithm toolbox
+ Heterogeneous positioning toolbox 
+ Various techniques exploiting heterogeneous Radio observables


Pylayers was developed over the time to help with research into impulse radio
propagation studies and  algorithms for
localization, but it may be also useful to help with the development of
characterization, Wireless Sensor Network, MANET, radio cooperation techniques, physical layer
security ...

The simulator is written in Python and C making use of NumPy, SciPy, SimPy,
Matplotlib, Scikit-Learn, CVXOPT, Shapely, Networkx. It is designed to be easy to
extend with new functionalities, allowing it to evolve over the time as a common
tool for research and development in a wide range of communications and
localization applications. We warmly encourage all users to contribute 
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

    manual/manual.rst 

    modules/index.rst
    
    bibliography.rst

Publications
------------

If you publish work which makes use of PyLayers, please cite the following paper:



Download and Installation
-------------------------

The current version is 0.1, released on 10 October 2011. Download the release in your preferred format:

Pylayers was developed on Linux and probably won't readily work on Windows platforms

There list of python dependencies is 




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

