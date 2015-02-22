from pylayers.antprop.coverage import *
C = Coverage()
C.cover()
f,a=C.show(typ='sinr',figsize=(10,8))
plt.show()
