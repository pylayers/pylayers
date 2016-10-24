# -*- coding: utf-8 -*-
r"""
=================================================
Building graphs of a Layout 
=================================================

"""
from pylayers.gis.layout import *
# Load the layout from its .ini file in $BASENAME/struc/ini 
L = Layout('WHERE1.ini')
# Build all the graphs 
# Check graphs
L._visual_check()
plt.show()

