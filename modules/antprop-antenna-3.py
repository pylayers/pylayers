import matplotlib.pyplot as plt
from pylayers.antprop.antenna import *
A = Antenna('S1R1.mat','ant/UWBAN/Matfile')
pol1 = plt.polar(A.phi,abs(A.Ftheta[10,45,:]),'b')
pol2 = plt.polar(A.phi,abs(A.Ftheta[20,45,:]),'r')
pol3 = plt.polar(A.phi,abs(A.Ftheta[30,45,:]),'g')
txt = plt.title('S1R1 antenna : st loadmat')
plt.show()
