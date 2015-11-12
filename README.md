# PyLayers 
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/logo.png" width="150">


PyLayers an Open Source Site Specific Radio Channel Simulator for 

    * Mobile Heterogeneous Networks,
    * Propagation,
    * Communications,
    * Mobility,
    * Localization,
    * Wearable Devices.

The code is developped at the University of Rennes 1 in the laboratory
Institute of Electronics and Telecommunication of Rennes. 

WebSite : [http://www.pylayers.org](http://www.pylayers.org)
Documentation : [Documentation](http://pylayers.github.io/pylayers/notebook/TOC.html)

Pylayers is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Main features

    * UWB Raytracing tool (Delays,DOA,DOD)
    * Indoor Radio Coverage
    * Human Mobility Simulator for Wearables and WBAN
    * Handling of Motion Capture
    * Rich Antenna Patterns Description
    * Heterogeneous Networks
    * Indoor Localization Platform
    * Handling of various Radio Standards including Ultra Wideband
    * ...

## Installation Notes

### Linux / OS X

The prefer way to install PyLayers is to first install the free python distribution [Anaconda](https://store.continuum.io/cshop/anaconda/) on your platform. 
By doing so, you're installing most of the required dependencies.

1. Download [Anaconda](https://store.continuum.io/cshop/anaconda/)
2. Install Anaconda with the given command ( $bash Anaconda-<version>-<plateform>.sh ). 
Important : Say yes to add Anaconda to your path
3. Clone PyLayers : **git clone https://github.com/pylayers/pylayers.git**
4. Run "./installer_unix" file from the PyLayers directory (YOUR FIRST NEED to add exe rights using the command: chmod +x ./installer_unix )
5. Done

### Windows Install


1. Download and Install [Anaconda](https://store.continuum.io/cshop/anaconda/) 
2. Download the Shapely package corresponding to your platform from 
[Unofficial Windows Binaries for Python Extension Packages #Shapely](http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely) *
3. Download the Basemap package corresponding to your platform from [Unofficial Windows Binaries for Python Extension Packages #Basemap](http://www.lfd.uci.edu/~gohlke/pythonlibs/#basemap)
4. Install Shapely and Basemap using the downloaded files
5. 3. Clone PyLayers : **git clone https://github.com/pylayers/pylayers.git**
6. Run installer.bat
7. Done


* Despite shapely is part of the Anaconda distribution, the Windows distribution doesn't contains the libgeos dependencies.

## First Run

```python
>>> from pylayers.simul.link import *
>>> DL = DLink()
>>> DL.eval()
>>> DL._show3()
```
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/_show3.png" width="300">

## Troubleshooting

If you have any trouble during the install process you should first refer to the INSTALL.TXT file.

For any problem during the use of the tool, please open an issue on the [github dedicated page](https://github.com/pylayers/pylayers/issues)

## Developpers

Pylayers is driven by professor [Bernard Uguen](mailto:bernard.uguen@univ-rennes1.fr) at [University of Rennes 1](www.univ-rennes1.fr), [IETR](www.ietr.fr) laboratory and [ESIR school of engineering](esir.univ-rennes1.fr)

Pylayers is currently developed at IETR by [Bernard Uguen](mailto:bernard.uguen@univ-rennes1.fr) , [Nicolas Amiot](mailto:nicolas.amiot@univ-rennes1.fr), Mohamed Laaraiedh and Meriem Mhedbhi, with the technical support of all the members from the Research Team of the Propagation and Localization team of the IETR (UMR CNRS 6164) lab.

Among early pylayers contributors : Friedman Tchoffo Talom, Louis Marie Aubert, Roxana Burghelea, Yu Lei, Taguhi Chaluyman, Stéphane Avrillon and Eric Plouhinec (Saint-Cyr CREC research center).

## Acknowledgements

This work has been supported by the Bretagne Region council (Project LOCUS), by the french ANR project AUBADE and CORMORAN and by the European projects FP7 UCELLS, FP7 WHERE1, FP7 WHERE2.

<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/logoUR1.jpg" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/logoIETR.jpg" width="100
">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/ucells.png" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/where1logo.jpg" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/WHERE2_Logo.jpg" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/Cormo.png" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/IR.png" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/bretagnegd.jpg" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/anr.png" width="100">
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/fp7.png" width="100">



License
Copyright ©, 2013 University of Rennes 1
<!-- 
Pylayers is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. -->






[![Analytics](https://ga-beacon.appspot.com/UA-34943220-2/pylayers/pylayers)](https://github.com/igrigorik/ga-beacon)
