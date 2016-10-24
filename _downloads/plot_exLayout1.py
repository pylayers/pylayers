# -*- coding: utf-8 -*-
r"""
=================================================
Loading an outdoor layout from its address
=================================================

"""
from pylayers.gis.layout import *
# Load the layout from open street map
L = Layout('Servon sur Vilaine',dist_m=50)
# Show the layout
L.showG('s')
plt.show()

