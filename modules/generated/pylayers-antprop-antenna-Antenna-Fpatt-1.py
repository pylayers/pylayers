from pylayers.antprop.antenna  import *
A = Antenna(typ='Gauss')
A.Fpatt()
f,a = A.polar()
plt.show()
