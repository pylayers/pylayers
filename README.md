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

See the [Install Page](https://github.com/pylayers/pylayers/blob/master/INSTALL.md)


## A Ray-Tracing simulation in 4 lines 

```python
>>> from pylayers.simul.link import *
>>> DL = DLink()
>>> DL.eval()
>>> DL._show3()
```
<img src="https://github.com/pylayers/pylayers/blob/master/doc/_static/_show3.png" width="300">

## Troubleshooting

If you have any trouble during the install process you should first refer to the [Install Page](https://github.com/pylayers/pylayers/blob/master/INSTALL.md) and particularly to the [Usual Install Issues](https://github.com/pylayers/pylayers/blob/master/INSTALL.md#usual-install-issues) part.

For any problem during the use of the tool, please open an issue on the [dedicated github's page](https://github.com/pylayers/pylayers/issues)

## Developpers

Pylayers is driven by professor [Bernard Uguen](mailto:bernard.uguen@univ-rennes1.fr) at [University of Rennes 1](www.univ-rennes1.fr), [IETR](www.ietr.fr) laboratory and [ESIR school of engineering](esir.univ-rennes1.fr)

Pylayers is currently developed at IETR by [Bernard Uguen](mailto:bernard.uguen@univ-rennes1.fr), [Nicolas Amiot](mailto:nicolas.amiot@univ-rennes1.fr) and Mamadou Dialounde Balde, with the technical support of all the members from the Research Team of the Propagation and Localization team of the IETR (UMR CNRS 6164) lab.

Among early pylayers contributors : Friedman Tchoffo Talom, Louis Marie Aubert, Roxana Burghelea, Yu Lei, Taguhi Chaluyman, Stéphane Avrillon, Mohamed Laaraiedh, Meriem Mhedbhi and Eric Plouhinec (Saint-Cyr CREC research center).

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
