import matplotlib.pyplot as plt
from pylayers.antprop.antenna import *
A = Antenna('mat','S1R1.mat','ant/UWBAN/Matfile')
pol = plt.polar(A.phi,abs(A.Ftheta[10,45,:]))
txt = plt.title('test loadmat')
plt.show()
