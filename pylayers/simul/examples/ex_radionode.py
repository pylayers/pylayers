from pylayers.simul.radionode import *

# create a RadioNode

tx = RadioNode()
tx.point([1,1,1],mode='append')
tx.point([1,2,1],mode='append')
tx.point([2,1,1],mode='append')
tx.info()

# visualize the RadioNode in geomview 

tx.show3()
