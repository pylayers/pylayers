pylayers
========

Pylayers is a simulator designed to model mobile pedestrian radio propagation
in Indoor environments.


The simulator is able to welcome various kind of radio propagation solvers from
the simpler to the more complex. The current version can generate UWB channel
impulse responses along trajectory thanks to external C functions which have
been developed few years ago in the french ANR project AUBADE.


Taking advantage of the project FP7 UCELLS, FP7 WHERE1 and FP7 WHERE2 the tool
has gradually been dedicated to localization applications. Pylayers contains
various functions which allows to estimate the position based a given number
of Location Dependent Parameters (LDPs)

Currently, the emphasis is

which are obtained by exploiting the Simpy dicrete event package and the
Transit code originaly developped by 

Features included in the current development version goes through different
layers.

+ Layout description and edition
+ Multi Layer transmission coefficient 
+ Account fot the full vector antenna pattern  
+ Motley Keenan model 
+ C coded Ray Tracing engine
+ Handling of Impulse Radio Ultra Wideband Waveforms 
+ Mobile user mobility 
+ Set Membership positoning algorithm toolbox
+ Heterogeneous positioning toolbox 
+ Various techniques exploiting heterogeneous Radio observables


Pylayers was developed over the time to help with research into algorithms for
localization, but may be also useful to help with the development of radion channel
characterization, Wireless Sensor Network, MANET, relaying, Physical Layer
Security.

The simulator is written in Python and C making use of NumPy, SciPy, SimPy,
Matplotlib, Scikit-Learn, Shapely, Networkx. It is designed to be easy to
extend with new functionality, allowing it to continuously evolve as a common
tool for research and development in a wide range of communications and
localization applications. We encourage all users to contribute new algorithms,
models and other improvements back to the project.

Just fork the github repository on github

IMUSim was written by Bernard Uguen, Nicolas Amiot, at the IETR, University of Rennes 1.
This work was supported by the Bretagne Region Council (Project LOCUS), by the french ANR
project AUBADE and CORMORAN and by project the european project FP7 UCELLS, FP7 WHERE1, FP7 WHERE2.


Publications
------------

If you publish work which makes use of PyLayers, please cite the following paper:



Download and Installation
-------------------------

The current version is 0.1, released on 10 October 2011. Download the release in your preferred format:

Pylayers was developed on Linux and probably won't work on Windows platforms

There is a large number of dependenceies and 
Documentation

<p>Please see the <a href="http://www.imusim.org/docs/index.html">documentation index for the current release</a>. There is a new <a href="http://www.imusim.org/docs/tutorial.html">tutorial</a> to help you get started with IMUSim.</p>

<h2>Mailing List</h2>

There is a "http://groups.google.com/group/pylayers-users" mailing list for users of Pylayers  on which you can get help with using or extending the software. New releases will also be announced on the mailing list. To join, enter your email address below.

http://groups.google.com/group/imusim-users/boxsubscribe">
Email: <input type="text" name="email"> <input type="submit" name="sub" value="Join">
</form><p></p>


License
-------

Copyright © 2012 University of Rennes 1

Pylayers is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Pylayers is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the <a href="http://www.gnu.org/licenses/">GNU General Public License for more details.

.. toctree::
    :maxdepth: 2 

    modules/index.rst
    
    UserManual.rst 

    bsignal.rst

    references.rst
