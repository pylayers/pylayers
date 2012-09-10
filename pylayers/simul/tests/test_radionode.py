from pylayers.simul.radionode import *

# create a RadioNode

tx = RadioNode(typ='tx')
tx.point([1,1,1])
tx.info()

# visualize the RadioNode in geomview 

tx.show3()
