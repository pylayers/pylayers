import matplotlib.pyplot as plt
from pylayers.antprop.antenna import *
A = Antenna('S1R1.mat',directory='ant/UWBAN/Matfile')
f,a = A.plotG(plan='theta',angdeg=0)
f,a = A.plotG(plan='phi',angdeg=90,fig=f,ax=a)
txt = plt.title('S1R1 antenna : st loadmat')
plt.show()
